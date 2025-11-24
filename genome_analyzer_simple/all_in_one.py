#!/usr/bin/env python3
"""
게놈 분석기 코어 (Genome Analyzer Core - v4.0)
통합: KNIH rpoB 공식 스크립트 연동 + MLST + AMR/Plasmid 탐지

[핵심 변경 사항]
1. 아키텍처 최적화 (Wrapper Pattern): rpoB 분석 시, 복잡한 재구현 대신 
   제공된 'identify_from_genome.py' 스크립트를 서브프로세스로 호출합니다.
2. 데이터 파이프라인: 외부 스크립트가 생성한 TSV 결과를 파싱하여 
   전체 통합 리포트(Final_Analysis_Report.txt)에 병합합니다.
"""

import asyncio
import argparse
import subprocess
import shutil
import json
import time
import sys
import os
from pathlib import Path
from typing import Tuple, List, Dict, Any
import pandas as pd

# --- [Level 0] 설정 상수 ---
DEFAULT_IDENTITY = 90.0
DEFAULT_COVERAGE = 90.0
REQUIRED_TOOLS = ["blastn", "makeblastdb", "blastdbcmd"]

# ==============================================================================
# [Level 1] 메인
# ==============================================================================

async def main():
    parser = argparse.ArgumentParser(description="Genome Analyzer Enterprise Core v4")
    parser.add_argument("-g", "--genome", required=True, type=Path, help="입력 FASTA 파일 경로")
    parser.add_argument("-s", "--species", required=True, help="타겟 종 (MLST DB 매칭)")
    parser.add_argument("-o", "--output", type=Path, default=Path("analysis_results"), help="결과물 저장 경로")
    parser.add_argument("--db", type=Path, default=Path("database"), help="참조 데이터베이스 루트 (예: ./db)")
    parser.add_argument("--audit-mode", action="store_true", help="검증을 위한 중간 분석 데이터(Artifacts) 보존")
    args = parser.parse_args()

    check_system_requirements(args.db)

    # 작업 초기화
    start_time = time.time()
    genome_id = args.genome.stem
    res_dir = args.output / genome_id
    temp_dir = res_dir / "audit_artifacts"
    
    res_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(exist_ok=True)
    
    log_info(f"분석 세션 시작: {genome_id}")

    try:
        # 1. 게놈 인덱싱 (Blocking) - MLST 및 Feature 탐지에 필요
        genome_db = await create_genome_db(args.genome, temp_dir)

        # 2. 작업 스케줄링 (Parallel Execution)
        tasks = []
        
        # [Task 1] rpoB 정밀 종 동정 (KNIH 스크립트 래퍼)
        # 외부 스크립트 탐색: --db 경로와 같은 레벨 혹은 상위에 있는지 확인
        rpoB_script = find_rpoB_script(args.db)
        
        if rpoB_script:
            # 외부 스크립트는 원본 FASTA 파일(args.genome)을 필요로 함
            tasks.append(run_rpoB_wrapper(rpoB_script, args.genome, args.db, res_dir))
        else:
            log_warn("'identify_from_genome.py'를 찾을 수 없어 rpoB 분석을 건너뜁니다.")
            tasks.append(asyncio.create_task(asyncio.sleep(0, result={})))

        # [Task 2] MLST (분자 역학)
        mlst_db_root = args.db.parent if args.db.name == "db" else args.db
        if (mlst_db_root / "MLST_DB").exists():
             tasks.append(run_mlst(args.species, genome_db, mlst_db_root, res_dir, temp_dir))
        else:
             tasks.append(run_mlst(args.species, genome_db, args.db, res_dir, temp_dir))
        
        # [Task 3, 4] AMR / Plasmid
        for task_name, folder_name in [("Antimicrobial_Resistance", "resfinder_db"), ("Plasmid_Replicons", "plasmidfinder_db")]:
            target_db = None
            if (args.db / folder_name).exists(): target_db = args.db / folder_name
            elif (args.db.parent / folder_name).exists(): target_db = args.db.parent / folder_name
            
            if target_db:
                tasks.append(run_feature_detection(task_name, target_db, genome_db, res_dir))

        # 3. 실행 및 결과 집계
        results = await asyncio.gather(*tasks)
        
        rpoB_res = results[0]
        mlst_res = results[1]
        
        features_res = []
        for res_list in results[2:]:
            if isinstance(res_list, list):
                features_res.extend(res_list)
        
        # 4. 리포트 작성
        generate_report(genome_id, res_dir, rpoB_res, mlst_res, features_res)

    except Exception as e:
        log_error(f"치명적 런타임 오류: {e}")
        import traceback
        traceback.print_exc()
    finally:
        finalize_cleanup(temp_dir, args.audit_mode)
            
    log_info(f"작업 완료 (소요 시간: {time.time() - start_time:.2f}초)")


# ==============================================================================
# [Level 2]
# ==============================================================================

async def run_rpoB_wrapper(script_path: Path, genome_fasta: Path, db_dir: Path, output_dir: Path) -> Dict:
    """
    [rpoB Wrapper 모듈]
    외부 스크립트(identify_from_genome.py)를 실행하고 결과를 파싱합니다.
    """
    log_info(f"Task 시작: rpoB 정밀 동정 (External Script)")
    rpoB_out = output_dir / "rpoB_Identification"
    rpoB_out.mkdir(parents=True, exist_ok=True)

    # 외부 스크립트 실행 커맨드 구성
    cmd = [
        sys.executable, str(script_path),
        "--genome", str(genome_fasta),
        "--db-dir", str(db_dir),
        "--outdir", str(rpoB_out),
        "--threads", "4" # 병렬 처리 활용
    ]

    success, stdout, stderr = await run_command(cmd)
    
    if not success:
        log_error(f"rpoB 스크립트 실행 실패: {stderr}")
        return {"status": "SCRIPT_FAILED", "error": stderr}

    # 결과 파일 파싱 (diversity_hits.tsv)
    # 외부 스크립트가 생성하는 파일명을 기준으로 함
    result_tsv = rpoB_out / "diversity_hits.tsv"
    top_hits = []
    
    if result_tsv.exists():
        try:
            # 외부 스크립트의 컬럼 순서: qseqid, sseqid, pident, length, ...
            cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
            df = pd.read_csv(result_tsv, sep='\t', names=cols)
            
            if not df.empty:
                # 상위 5개 추출
                for _, row in df.head(5).iterrows():
                    top_hits.append({
                        "Reference_ID": row['sseqid'],
                        "Identity_Pct": row['pident'],
                        "Coverage_Pct": "N/A" # 외부 스크립트 출력에 따라 조정 필요 (여기서는 qlen/length 등으로 추정 가능하나 생략)
                    })
        except Exception as e:
            log_warn(f"rpoB 결과 파싱 중 오류: {e}")

    return {
        "status": "SUCCESS",
        "top_hit": top_hits[0]['Reference_ID'] if top_hits else "Unknown",
        "details": top_hits
    }

def find_rpoB_script(db_path: Path) -> Path:
    """
    [Utils] identify_from_genome.py 스크립트 위치 탐색
    DB 폴더와 같은 레벨, 혹은 상위 레벨, 혹은 scripts 폴더 내부 등을 탐색
    """
    candidates = [
        db_path / "identify_from_genome.py",
        db_path.parent / "identify_from_genome.py",
        db_path.parent / "scripts" / "identify_from_genome.py",
        Path("./identify_from_genome.py"),
        Path("./scripts/identify_from_genome.py")
    ]
    
    for path in candidates:
        if path.exists():
            return path.resolve()
    return None

async def run_mlst(species: str, genome_db: Path, db_root: Path, output_dir: Path, temp_dir: Path) -> Dict:
    """[MLST 모듈] 기존 로직 유지"""
    log_info(f"Task 시작: MLST ({species})")
    mlst_out = output_dir / "MLST"
    mlst_out.mkdir(parents=True, exist_ok=True)
    
    species_db = db_root / "MLST_DB" / species
    profile_file = species_db / f"{species}.txt"
    
    if not profile_file.exists():
        log_error(f"프로파일 정의 파일 누락: {profile_file}")
        return {"status": "PROFILE_ERROR"}

    with open(profile_file, 'r') as f:
        header = f.readline().strip().split('\t')
        loci = header[1:]

    probes = temp_dir / "loci_probes.fasta"
    with open(probes, 'w') as f_out:
        for locus in loci:
            locus_file = species_db / f"{locus}.tfa"
            if locus_file.exists():
                f_out.write(locus_file.read_text())
    
    probe_hits = temp_dir / "loci_mapping.tsv"
    await run_blastn(query=probes, db=genome_db, out=probe_hits, options={"perc_identity": 90})
    
    extracted_genes = temp_dir / "target_loci_extracted.fasta"
    try:
        df = pd.read_csv(probe_hits, sep='\t', names=['q', 's', 'pident', 'len', 'mis', 'gap', 'qs', 'qe', 'ss', 'se', 'e', 'bit', 'cov'])
        best_hits = df.loc[df.groupby('q')['bitscore'].idxmax()] 
    except pd.errors.EmptyDataError:
        return {"st": "DETECTION_FAILED"}

    processed_loci = set()
    with open(extracted_genes, 'w') as f_ext:
        for _, row in best_hits.iterrows():
            locus_name = row['q'].split('_')[0]
            if locus_name in processed_loci: continue
            s, e = sorted((int(row['ss']), int(row['se'])))
            strand = "plus" if int(row['ss']) < int(row['se']) else "minus"
            cmd = ["blastdbcmd", "-db", str(genome_db), "-entry", row['s'], "-range", f"{s}-{e}", "-strand", strand]
            succ, seq_data, _ = await run_command(cmd)
            if succ:
                seq_lines = seq_data.strip().split('\n')
                if len(seq_lines) > 1:
                    f_ext.write(f">{locus_name}\n{''.join(seq_lines[1:])}\n")
                    processed_loci.add(locus_name)

    all_alleles = temp_dir / "allele_library.fasta"
    with open(all_alleles, 'w') as f_out:
        for locus_file in species_db.glob("*.tfa"):
            f_out.write(locus_file.read_text())
            
    allele_db = await create_genome_db(all_alleles, temp_dir)
    allele_calls = temp_dir / "allele_identification.tsv"
    await run_blastn(query=extracted_genes, db=allele_db, out=allele_calls, options={"perc_identity": 99, "outfmt": "6 qseqid sseqid pident length"})
    
    profile_map = {}
    try:
        df_a = pd.read_csv(allele_calls, sep='\t', names=['q', 's', 'pident', 'len'])
        df_a = df_a.sort_values('pident', ascending=False)
        for locus in loci:
            hit = df_a[df_a['q'] == locus]
            if not hit.empty:
                top = hit.iloc[0]
                if top['pident'] == 100.0:
                    import re
                    match = re.search(r'_(\d+)$', top['s'])
                    profile_map[locus] = match.group(1) if match else "?"
                else:
                    profile_map[locus] = f"~{top['s']}({top['pident']}%)" 
            else:
                profile_map[locus] = "-"
    except:
        return {"st": "PROCESSING_ERROR"}

    df_profile = pd.read_csv(profile_file, sep='\t').astype(str)
    query_profile = [profile_map.get(l, "-") for l in loci]
    found_st = "Unknown (Novel)"
    for _, row in df_profile.iterrows():
        if query_profile == [str(row[l]) for l in loci]:
            found_st = row['ST']
            break
            
    result = {"st": found_st, "profile": profile_map}
    with open(mlst_out / "mlst_report.json", 'w') as f:
        json.dump(result, f, indent=4)
    return result

async def run_feature_detection(analysis_name: str, db_folder: Path, genome_db: Path, output_dir: Path) -> List[Dict]:
    """[특징 탐지 모듈] AMR, Plasmid 등 일반 참조 검색"""
    log_info(f"Task 시작: {analysis_name}")
    analysis_dir = output_dir / analysis_name
    analysis_dir.mkdir(parents=True, exist_ok=True)
    
    refs = list(db_folder.rglob("*.f*a"))
    if not refs: return []

    combined_ref = analysis_dir / "reference_library.fasta"
    with open(combined_ref, 'w') as f_out:
        for ref in refs: f_out.write(ref.read_text())

    blast_out = analysis_dir / "detection_results.tsv"
    await run_blastn(query=combined_ref, db=genome_db, out=blast_out, options={"perc_identity": DEFAULT_IDENTITY, "qcov_hsp_perc": DEFAULT_COVERAGE})

    results = []
    try:
        cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qs', 'qe', 'ss', 'se', 'e', 'bit', 'cov']
        df = pd.read_csv(blast_out, sep='\t', names=cols)
        if not df.empty:
            best_hits = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]
            for _, row in best_hits.iterrows():
                results.append({
                    "Type": analysis_name,
                    "Gene": row['qseqid'],
                    "Contig_ID": row['sseqid'],
                    "Identity_Pct": row['pident'],
                    "Coverage_Pct": row['cov']
                })
    except: pass
    with open(analysis_dir / f"{analysis_name}_summary.json", 'w') as f:
        json.dump(results, f, indent=4)
    return results

# ==============================================================================
# [Level 3]
# ==============================================================================

async def create_genome_db(fasta_path: Path, output_dir: Path) -> Path:
    """[Engine] makeblastdb 래퍼 (캐싱 지원)"""
    db_name = output_dir / fasta_path.stem
    if not all(db_name.with_suffix(s).exists() for s in ['.nin', '.nhr', '.nsq']):
        log_info(f"DB 인덱싱 수행: {fasta_path.name}")
        cmd = ["makeblastdb", "-in", str(fasta_path), "-dbtype", "nucl", "-out", str(db_name), "-parse_seqids"]
        success, _, stderr = await run_command(cmd)
        if not success: raise RuntimeError(f"makeblastdb 실패: {stderr}")
    return db_name

async def run_blastn(query: Path, db: Path, out: Path, options: Dict[str, Any] = None):
    """[Engine] blastn 래퍼"""
    opts = {"outfmt": "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp"}
    if options: opts.update(options)
    cmd = ["blastn", "-query", str(query), "-db", str(db), "-out", str(out)]
    for k, v in opts.items(): cmd.extend([f"-{k}", str(v)])
    success, _, stderr = await run_command(cmd)
    if not success: raise RuntimeError(f"blastn 실패: {stderr}")

def generate_report(genome_id: str, results_dir: Path, rpoB_data: Dict, mlst_data: Dict, features_data: List):
    """[Utils] 통합 리포트 생성기"""
    report_file = results_dir / "Final_Analysis_Report.txt"
    with open(report_file, 'w') as f:
        f.write("=== 게놈 분석 최종 리포트 (Genome Analysis Final Report) ===\n")
        f.write(f"Sample ID: {genome_id}\n")
        f.write(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("--- 1. rpoB 종 동정 (Species Identification) ---\n")
        if rpoB_data and rpoB_data.get("status") == "SUCCESS":
            f.write(f"최우선 동정 결과: {rpoB_data.get('top_hit')}\n")
            f.write("상위 후보 리스트 (Top 5 Candidates):\n")
            for hit in rpoB_data.get('details', []):
                f.write(f"  - {hit['Reference_ID']} (Identity: {hit['Identity_Pct']}%)\n")
        else:
            f.write(f"동정 실패 또는 데이터 없음. 상태: {rpoB_data.get('status', 'N/A')}\n")
            if rpoB_data.get("error"):
                f.write(f"Error Log: {rpoB_data.get('error')}\n")
        f.write("\n")

        f.write("--- 2. 분자 역학 분석 (MLST) ---\n")
        f.write(f"Sequence Type (ST): {mlst_data.get('st', 'N/A')}\n")
        f.write(f"Allele Profile: {json.dumps(mlst_data.get('profile', {}))}\n\n")
        
        f.write("--- 3. 유전자 마커 탐지 결과 ---\n")
        if not features_data:
            f.write("검출된 마커 없음.\n")
        else:
            f.write(f"{'Type':<20} {'Gene':<20} {'Contig':<20} {'Identity(%)':<12} {'Coverage(%)'}\n")
            f.write("-" * 85 + "\n")
            for hit in features_data:
                f.write(f"{hit['Type']:<20} {hit['Gene']:<20} {hit['Contig_ID']:<20} {hit['Identity_Pct']:<12} {hit['Coverage_Pct']}\n")
    
    log_info(f"리포트 생성 완료: {report_file}")

def check_system_requirements(db_path: Path):
    for tool in REQUIRED_TOOLS:
        if not shutil.which(tool):
            log_error(f"시스템 구성 오류: 필수 도구 '{tool}' 미설치")
            sys.exit(1)
    if not db_path.exists():
        log_error(f"데이터베이스 로드 실패: 경로 '{db_path}' 없음")
        sys.exit(1)

def finalize_cleanup(temp_dir: Path, audit_mode: bool):
    if not audit_mode and temp_dir.exists():
        shutil.rmtree(temp_dir)
        log_info("임시 아티팩트 정리 완료.")
    elif audit_mode:
        log_info(f"감사 추적 데이터 보존됨: {temp_dir}")

# ==============================================================================
# [Level 4] 시스템 I/O
# ==============================================================================

async def run_command(cmd: List[str]) -> Tuple[bool, str, str]:
    cmd_str = [str(c) for c in cmd]
    try:
        proc = await asyncio.create_subprocess_exec(*cmd_str, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE)
        stdout, stderr = await proc.communicate()
        return proc.returncode == 0, stdout.decode('utf-8', errors='ignore'), stderr.decode('utf-8', errors='ignore')
    except Exception as e: return False, "", str(e)

def log_info(msg: str): print(f"[INFO] {time.strftime('%H:%M:%S')} | {msg}")
def log_warn(msg: str): print(f"[WARN] {time.strftime('%H:%M:%S')} | {msg}")
def log_error(msg: str): print(f"[ERROR] {time.strftime('%H:%M:%S')} | {msg}", file=sys.stderr)

if __name__ == "__main__":
    asyncio.run(main())