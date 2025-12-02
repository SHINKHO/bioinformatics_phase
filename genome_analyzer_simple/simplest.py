import argparse
import asyncio
import json
import os
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd

# --- 시스템 상수 ---
# 필수 외부 도구 목록
REQUIRED_TOOLS = ["blastn", "makeblastdb", "blastdbcmd"]
# BLAST 실행 시 사용할 스레드 수
NUM_THREADS = "4"
# BLAST 검색 시 기본 Identity (유사도) 임계값
DEFAULT_IDENTITY = 90.0
# BLAST 검색 시 기본 Coverage (커버리지) 임계값
DEFAULT_COVERAGE = 90.0


# ==============================================================================
# 메인 실행 함수
# ==============================================================================

async def main():
    """메인 함수: 전체 유전체 분석 파이프라인을 실행합니다."""
    # 1. 인자 파싱: 커맨드라인에서 입력받은 인자를 처리합니다.
    parser = argparse.ArgumentParser(description="Genome Analyzer Direct Core")
    parser.add_argument("-g", "--genome", required=True, type=Path, help="Input FASTA file")
    parser.add_argument("-s", "--species", required=True, help="Target species for MLST")
    parser.add_argument("-o", "--output", type=Path, default=Path("analysis_results"), help="Output directory")
    parser.add_argument("--db", type=Path, default=Path("database"), help="Database root directory")
    parser.add_argument("--audit-mode", action="store_true", help="debugging : 중간 산출물 보존")
    args = parser.parse_args()

    # 2. 필수 도구 확인
    _check_required_tools()

    # 3. 환경 설정: 시간 측정 시작 및 결과 저장 경로 설정
    start_time = time.time()
    genome_id = args.genome.stem
    res_dir = args.output / genome_id
    temp_dir = res_dir / "artifacts"

    res_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(exist_ok=True)

    print(f"[SYSTEM] Analysis session started: {genome_id} (Threads: {NUM_THREADS})")

    try:
        # 4. 게놈 DB 생성: 입력된 유전체로 BLAST DB를 만듭니다.
        genome_db_prefix = temp_dir / genome_id
        await make_blast_db(args.genome, genome_db_prefix)

        # 5. 균종동정 : rpoB 처리하여 결과값을 받아 분기 처리합니다
        rpoB_result = {}  
        rpoB_script = find_script_smart(args.db, "identify_from_genome.py", "rpoB_db")
        
        if rpoB_script:
            rpoB_db_dir = rpoB_script.parent / "db"
            if not rpoB_db_dir.exists():
                rpoB_db_dir = args.db
            rpoB_result = await task_rpoB_wrapper(rpoB_script, args.genome, rpoB_db_dir, res_dir)
        else:
            print("[ERROR] rpoB script not found. failed species identification.")
            sys.exit(0)
        
        # rpoB 결과가 Unknown이면 종료
        if rpoB_result.get("species", "Unknown") == "Unknown":
            print("[WARN] rpoB Process. failed species identification. Aborting.")
            sys.exit(0)

        # 5. 분석 작업 생성: 모든 분석을 비동기 태스크로 만듭니다.
        tasks = _create_analysis_tasks(args, genome_db_prefix, res_dir, temp_dir)

        # 6. 모든 분석 동시 실행
        results = await asyncio.gather(*tasks)

        # 7. 결과 취합
        data_collection = {
            "genome_id": genome_id,
            "species_rpoB": rpoB_result,
            "mlst": results[0],
            "amr": results[1],
            "plasmid": results[2],
            "virulence": results[3],
            "mge": results[4],
        }

        # 8. 최종 보고서 작성
        write_structured_report(res_dir / "Final_Report.txt", data_collection)

    except Exception as e:
        # 런타임 오류 처리
        print(f"[ERROR] A runtime error occurred: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
    finally:
        # 임시 파일 정리
        if not args.audit_mode and temp_dir.exists():
            shutil.rmtree(temp_dir)
            print("[SYSTEM] Temporary files cleaned up.")

    # 전체 실행 시간 출력
    print(f"[SYSTEM] Completed in {time.time() - start_time:.2f} seconds.")


def _check_required_tools():
    """필수 커맨드라인 도구들이 시스템에 설치되어 있는지 확인합니다."""
    for tool in REQUIRED_TOOLS:
        if not shutil.which(tool):
            print(f"[CRITICAL] Required tool not found: {tool}", file=sys.stderr)
            sys.exit(1)


def _create_analysis_tasks(args: argparse.Namespace, genome_db_prefix: Path, res_dir: Path, temp_dir: Path) -> List[
    asyncio.Task]:
    """모든 분석에 대한 비동기 태스크 리스트를 생성합니다."""
    tasks = []

    # 작업 1: MLST
    tasks.append(task_mlst(args.species, genome_db_prefix, args.db, res_dir, temp_dir))

    # 작업 2 ~: 유전자 특성 검색 (AMR, Plasmid 등)
    feature_searches = [
        ("Antimicrobial_Resistance", "resfinder_db"),
        ("Plasmid_Replicons", "plasmidfinder_db"),
        ("Virulence_Factors", "virulencefinder_db", "Pathogenfinder"),
        ("Mobile_Genetic_Elements", "MEFinder"),
    ]

    for search in feature_searches:
        name, db_name, *fallback = search
        db_path = args.db / db_name
        if not db_path.exists() and fallback:
            db_path = args.db / fallback[0]

        if db_path.exists():
            tasks.append(task_feature_search(name, db_path, genome_db_prefix, res_dir))
        else:
            print(f"[WARN] {name} database not found at '{db_path}'. Skipping.")
            tasks.append(asyncio.create_task(asyncio.sleep(0, result=[])))

    return tasks


# ==============================================================================
# 분석 태스크 래퍼 함수
# ==============================================================================

async def make_blast_db(fasta_input: Path, db_prefix: Path):
    """FASTA 파일로부터 BLAST 데이터베이스를 생성합니다 (이미 존재하지 않을 경우)."""
    # DB 파일들이 이미 존재하면 함수 종료
    if all(db_prefix.with_suffix(s).exists() for s in ['.nin', '.nhr', '.nsq']):
        return
    # makeblastdb 명령어 실행
    cmd = ["makeblastdb", "-in", str(fasta_input), "-dbtype", "nucl", "-out", str(db_prefix), "-parse_seqids"]
    success, _, stderr = await execute_command(cmd)
    if not success:
        raise RuntimeError(f"makeblastdb failed for {fasta_input}: {stderr}")


async def task_rpoB_wrapper(script_path: Path, genome_fa: Path, db_dir: Path, out_dir: Path) -> Dict[str, Any]:
    """rpoB 종 동정 스크립트를 실행하는 래퍼 함수입니다."""
    print("[TASK] Starting rpoB analysis")

    # 결과 저장 경로 생성
    final_out_dir = out_dir / "rpoB_Identification"
    final_out_dir.mkdir(exist_ok=True, parents=True)

    # 스크립트 실행 환경 변수 설정
    env = os.environ.copy()
    env["RPOB_DB_DIR"] = str(db_dir.resolve())
    ref_core_path = db_dir.parent / "ref" / "rpoB_core.fna"
    if ref_core_path.exists():
        env["RPOB_REF_CORE"] = str(ref_core_path.resolve())

    # 외부 스크립트 실행
    cmd = [sys.executable, str(script_path.resolve()), str(genome_fa.resolve())]
    proc = await asyncio.create_subprocess_exec(
        *cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE,
        cwd=script_path.parent, env=env
    )
    _, stderr = await proc.communicate()

    if proc.returncode != 0:
        return {"status": "FAILED", "error": stderr.decode('utf-8', errors='ignore')}

    # 결과 파일 이동 및 최상위 히트 파싱
    top_hit = "Unknown"
    generated_results_dir = script_path.parent / "results"
    if generated_results_dir.exists():
        for f in generated_results_dir.glob("*"):
            shutil.move(str(f), str(final_out_dir / f.name))
        try:
            generated_results_dir.rmdir()
        except OSError:
            pass

        tsv_path = final_out_dir / "diversity_hits.tsv"
        if tsv_path.exists() and tsv_path.stat().st_size > 0:
            try:
                df = pd.read_csv(tsv_path, sep='\t', header=None)
                if not df.empty and df.shape[1] >= 2:
                    raw_hit = str(df.iloc[0, 1])
                    top_hit = raw_hit.split("|")[0] if "|" in raw_hit else raw_hit
            except pd.errors.ParserError:
                pass

    return {"status": "SUCCESS", "species": top_hit}


async def task_mlst(species: str, genome_db: Path, db_root: Path, out_dir: Path, temp_dir: Path) -> Dict[str, Any]:
    """MLST 분석을 수행합니다."""
    print(f"[TASK] Starting MLST analysis for {species}")
    species_dir = db_root / "MLST_DB" / species

    # 1. 프로파일 및 유효한 Loci(유전자 좌위) 찾기
    profile_path = species_dir / f"{species}.txt"
    if not profile_path.exists():
        found = list(species_dir.glob("*.txt"))
        if not found:
            return {"st": "Unknown (Profile Missing)", "profile": {}}
        profile_path = found[0]

    with open(profile_path, 'r') as f:
        raw_headers = f.readline().strip().split('\t')
    potential_loci = raw_headers[1:]
    valid_loci = [l for l in potential_loci if (species_dir / f"{l}.tfa").exists()]

    if not valid_loci:
        return {"st": "Unknown (No Valid Loci Found)", "profile": {}}
    print(f"       -> Valid Loci: {valid_loci}")

    # 2. 모든 Loci에 대한 프로브(probe) 서열을 하나의 파일로 합치기
    probes_fa = temp_dir / "mlst_probes.fasta"
    with open(probes_fa, 'w') as f_out:
        for locus in valid_loci:
            f_out.write((species_dir / f"{locus}.tfa").read_text())

    # 3. BLAST (1차): 게놈 DB에 프로브 서열을 매핑
    map_tsv = temp_dir / "mlst_map.tsv"
    blast_cmd_map = [
        "blastn", "-query", str(probes_fa), "-db", str(genome_db), "-out", str(map_tsv),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-perc_identity", "95", "-num_threads", NUM_THREADS
    ]
    success, _, stderr = await execute_command(blast_cmd_map)
    if not success:
        print(f"[ERROR] MLST blastn (mapping) failed: {stderr}", file=sys.stderr)
        return {"st": "Unknown (Mapping Fail)", "profile": {}}

    # 4. BLAST 결과 파싱 및 최적 히트 추출
    try:
        df = pd.read_csv(map_tsv, sep='\t',
                         names=['q', 's', 'id', 'len', 'mis', 'gap', 'qs', 'qe', 'ss', 'se', 'e', 'bit'])
        if df.empty:
            return {"st": "Unknown (No Hits)", "profile": {}}
        best_hits = df.loc[df.groupby('q')['bit'].idxmax()]
    except pd.errors.EmptyDataError:
        return {"st": "Unknown (No Hits)", "profile": {}}
    except Exception as e:
        print(f"[ERROR] MLST mapping result processing failed: {e}", file=sys.stderr)
        return {"st": "Unknown (Detection Fail)", "profile": {}}

    # 5. 게놈에서 매칭된 서열 추출
    extracted_fa = temp_dir / "mlst_extracted.fasta"
    processed = set()
    with open(extracted_fa, 'w') as f_ext:
        for _, row in best_hits.iterrows():
            locus = row['q'].split('_')[0]
            if locus in processed:
                continue

            start, end = sorted((int(row['ss']), int(row['se'])))
            strand = "plus" if int(row['ss']) < int(row['se']) else "minus"

            cmd_extract = ["blastdbcmd", "-db", str(genome_db), "-entry", row['s'], "-range", f"{start}-{end}",
                           "-strand", strand]
            succ, seq, _ = await execute_command(cmd_extract)

            if succ:
                seq_cleaned = ''.join(seq.strip().splitlines()[1:])
                f_ext.write(f">{locus}\n{seq_cleaned}\n")
                processed.add(locus)

    # 6. BLAST (2차): 추출된 서열을 모든 대립유전자(allele) DB와 비교하여 동정
    all_alleles_fa = temp_dir / "mlst_all_alleles.fasta"
    with open(all_alleles_fa, 'w') as f_out:
        for locus in valid_loci:
            f_out.write((species_dir / f"{locus}.tfa").read_text())

    allele_db_idx = temp_dir / "mlst_allele_db"
    await make_blast_db(all_alleles_fa, allele_db_idx)

    ident_tsv = temp_dir / "mlst_ident.tsv"
    blast_cmd_ident = [
        "blastn", "-query", str(extracted_fa), "-db", str(allele_db_idx), "-out", str(ident_tsv),
        "-outfmt", "6 qseqid sseqid pident length", "-perc_identity", "100", "-num_threads", NUM_THREADS
    ]
    success, _, stderr = await execute_command(blast_cmd_ident)
    if not success:
        print(f"[ERROR] MLST blastn (identification) failed: {stderr}", file=sys.stderr)

    # 7. 동정 결과로 프로파일 맵 생성
    profile_map = {}
    try:
        df_id = pd.read_csv(ident_tsv, sep='\t', names=['q', 's', 'id', 'len'])
        df_id = df_id.sort_values('id', ascending=False)
        for locus in valid_loci:
            hit = df_id[df_id['q'] == locus]
            if not hit.empty:
                top_hit = hit.iloc[0]
                match = re.search(r'_(\d+)$', top_hit['s'])
                allele_num = match.group(1) if match else "?"
                profile_map[locus] = allele_num if top_hit['id'] == 100.0 else f"~{allele_num}"
            else:
                profile_map[locus] = "-"
    except pd.errors.EmptyDataError:
        pass

    # 8. 프로파일 맵을 이용해 최종 ST(Sequence Type) 결정
    st = "Unknown"
    try:
        df_prof = pd.read_csv(profile_path, sep='\t').astype(str)
        query_vec = [profile_map.get(l, "-") for l in valid_loci]

        for _, row in df_prof.iterrows():
            db_vec = [str(row[l]) for l in valid_loci]
            if query_vec == db_vec:
                st = row.iloc[0]
                break
    except Exception as e:
        print(f"[ERROR] MLST Profile Matching Failed: {e}")

    return {"st": st, "profile": profile_map}


async def task_feature_search(name: str, db_path: Path, genome_db: Path, out_dir: Path) -> List[Dict[str, Any]]:
    """BLAST를 사용하여 유전자, 플라스미드 등 특정 유전적 특성을 찾는 범용 태스크입니다."""
    print(f"[TASK] Starting {name} search")
    task_dir = out_dir / name
    task_dir.mkdir(exist_ok=True)

    # 1. DB 경로에서 참조 서열 파일 수집
    refs = list(db_path.rglob("*.[fF][sSna][aA]"))
    if not refs:
        refs = [p for p in db_path.glob("*") if p.is_file() and p.suffix not in ['.txt', '.py', '.sh', '.md']]

    # "all."로 시작하는 파일은 제외 (중복 방지)
    filtered_refs = [r for r in refs if not r.name.startswith("all.")]
    if not filtered_refs:
        print(f"[WARN] No database files found for {name}.")
        return []

    # 2. 모든 참조 서열을 하나의 파일로 합치기
    combined_fa = task_dir / "combined_ref.fasta"
    with open(combined_fa, 'w') as f_out:
        for r_path in filtered_refs:
            content = r_path.read_text()
            f_out.write(content)
            if not content.endswith("\n"):
                f_out.write("\n")

    # 3. BLAST 실행: 합쳐진 참조 서열을 게놈 DB에 검색
    res_tsv = task_dir / "results.tsv"
    blast_cmd = [
        "blastn", "-query", str(combined_fa), "-db", str(genome_db), "-out", str(res_tsv),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp",
        "-perc_identity", str(DEFAULT_IDENTITY), "-qcov_hsp_perc", str(DEFAULT_COVERAGE),
        "-num_threads", NUM_THREADS,
    ]
    success, _, stderr = await execute_command(blast_cmd)
    if not success:
        print(f"[ERROR] {name} blastn failed: {stderr}", file=sys.stderr)
        return []

    # 4. BLAST 결과 처리
    hits = []
    try:
        cols = ['q', 's', 'id', 'len', 'mis', 'gap', 'qs', 'qe', 'ss', 'se', 'e', 'bit', 'cov']
        df = pd.read_csv(res_tsv, sep='\t', names=cols)
        if not df.empty:
            # 유전자 이름 정규화 (예: "tet(A)_1_AB123" -> "tet(A)")
            df['GeneName'] = df['q'].apply(lambda x: x.split('_')[0])
            # 동일 유전자에 대해 가장 점수가 높은 히트만 선택 (중복 제거)
            best_hits = df.sort_values(['GeneName', 'bit', 'id'], ascending=[True, False, False]).drop_duplicates(
                'GeneName', keep='first')

            # 최종 결과 리스트에 추가
            for _, row in best_hits.iterrows():
                hits.append({
                    "Gene": row['GeneName'],
                    "FullHeader": row['q'],
                    "Identity": row['id'],
                    "Coverage": row['cov'],
                })
    except pd.errors.EmptyDataError:
        pass

    # 5. 결과를 JSON 파일로 저장
    with open(task_dir / f"{name}.json", 'w') as f_json:
        json.dump(hits, f_json, indent=4)

    return hits


# ==============================================================================
# 유틸리티 함수
# ==============================================================================

async def execute_command(cmd: List[str], cwd: Path = None) -> Tuple[bool, str, str]:
    """외부 커맨드를 실행하고 (성공 여부, stdout, stderr)를 반환합니다."""
    cmd_str = [str(c) for c in cmd]
    try:
        proc = await asyncio.create_subprocess_exec(
            *cmd_str, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE, cwd=cwd
        )
        out, err = await proc.communicate()
        return proc.returncode == 0, out.decode('utf-8', errors='ignore'), err.decode('utf-8', errors='ignore')
    except Exception as e:
        return False, "", str(e)


def find_script_smart(base_path: Path, script_name: str, sub_dir: str = None) -> Path:
    """일반적인 위치에서 스크립트 파일을 찾습니다."""
    if sub_dir:
        if (base_path / sub_dir / script_name).exists():
            return (base_path / sub_dir / script_name).resolve()

    if (base_path / script_name).exists():
        return (base_path / script_name).resolve()

    found = list(base_path.rglob(script_name))
    return found[0].resolve() if found else None


def write_structured_report(path: Path, data: Dict[str, Any]):
    """최종 분석 결과를 구조화된 텍스트 파일로 작성합니다."""

    def format_hits(hits: List[Dict]) -> List[str]:
        """히트 리스트를 보고서 형식의 문자열 리스트로 변환합니다."""
        return [f"{h['Gene']} ({h['Identity']:.2f}%)" for h in hits] if hits else ["None"]

    with open(path, 'w') as f:
        f.write("== Analysis Report ==\n\n")

        f.write(f"Genome ID: {data['genome_id']}\n")

        species = data.get('species_rpoB', {}).get('species', 'Unknown')
        f.write(f"Species (rpoB-based): {species}\n\n")

        f.write("-- Molecular Epidemiology --\n")
        mlst_info = data.get('mlst', {})
        st = mlst_info.get('st', 'Unknown')
        f.write(f"  MLST: {st}\n")
        profile = mlst_info.get('profile', {})
        if profile:
            prof_str = ", ".join([f"{k}({v})" for k, v in profile.items()])
            f.write(f"  Allele Profile: {prof_str}\n\n")
        else:
            f.write("\n")

        f.write("-- Antimicrobial Resistance --\n")
        amr_genes = format_hits(data.get('amr', []))
        f.write("  Acquired Resistance Genes:\n")
        for gene in amr_genes:
            f.write(f"    - {gene}\n")


        f.write("-- Virulence Factors --\n")
        vir_genes = format_hits(data.get('virulence', []))
        for gene in vir_genes:
            f.write(f"  - {gene}\n")
        if not vir_genes or vir_genes[0] == "None":
            f.write("  Not Detected\n")
        f.write("\n")

        f.write("-- Mobile Genetic Elements --\n")
        plasmids = format_hits(data.get('plasmid', []))
        f.write("  Plasmids:\n")
        for plasmid in plasmids:
            f.write(f"    - {plasmid}\n")

        mges = format_hits(data.get('mge', []))
        f.write("  Other Mobile Genetic Elements:\n")
        for mge in mges:
            f.write(f"    - {mge}\n")


if __name__ == "__main__":
    # 스크립트가 직접 실행될 때 main 함수를 호출합니다.
    asyncio.run(main())