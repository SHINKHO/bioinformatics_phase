import asyncio
import argparse
import shutil
import json
import time
import sys
import os
import pandas as pd
import re
from pathlib import Path
from typing import List, Dict, Any

# --- 시스템 상수 ---
REQUIRED_TOOLS = ["blastn", "makeblastdb", "blastdbcmd"]
NUM_THREADS = "4"
DEFAULT_IDENTITY = 90.0
DEFAULT_COVERAGE = 90.0

# ==============================================================================
# 메인
# ==============================================================================

async def main():
    parser = argparse.ArgumentParser(description="Genome Analyzer Direct Core")
    parser.add_argument("-g", "--genome", required=True, type=Path, help="입력 FASTA 파일")
    parser.add_argument("-s", "--species", required=True, help="타겟 종 (MLST용)")
    parser.add_argument("-o", "--output", type=Path, default=Path("analysis_results"), help="결과 저장 경로")
    parser.add_argument("--db", type=Path, default=Path("database"), help="DB 루트 경로")
    parser.add_argument("--audit-mode", action="store_true", help="중간 산출물 보존 ")
    args = parser.parse_args()

    for tool in REQUIRED_TOOLS:
        if not shutil.which(tool):
            print(f"[CRITICAL] 필수 도구 누락: {tool}", file=sys.stderr)
            sys.exit(1)

    start_time = time.time()
    genome_id = args.genome.stem
    res_dir = args.output / genome_id
    temp_dir = res_dir / "artifacts"
    
    res_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(exist_ok=True)
    
    print(f"[SYSTEM] 분석 세션 시작: {genome_id} (Threads: {NUM_THREADS})")

    try:
        # 게놈 DB 인덱싱
        genome_db_prefix = temp_dir / genome_id
        await make_blast_db(args.genome, genome_db_prefix)

        # 분석 작업 리스트
        tasks = []

        # 1. rpoB (종 동정)
        fixed_rpoB_path = args.db / "rpoB_db" / "identify_from_genome.py"
        rpoB_script = fixed_rpoB_path if fixed_rpoB_path.exists() else find_script_smart(args.db, "identify_from_genome.py")
        
        if rpoB_script:
            rpoB_db_dir = rpoB_script.parent / "db"
            if not rpoB_db_dir.exists(): rpoB_db_dir = args.db
            tasks.append(task_rpoB_wrapper(rpoB_script, args.genome, rpoB_db_dir, res_dir))
        else:
            print("[WARN] rpoB 스크립트를 찾을 수 없습니다.")
            tasks.append(asyncio.create_task(asyncio.sleep(0, result={})))

        # 2. MLST
        tasks.append(task_mlst(args.species, genome_db_prefix, args.db, res_dir, temp_dir))

        # 3. AMR (ResFinder)
        amr_db = args.db / "resfinder_db"
        if amr_db.exists():
            tasks.append(task_feature_search("Antimicrobial_Resistance", amr_db, genome_db_prefix, res_dir))
        else:
            tasks.append(asyncio.create_task(asyncio.sleep(0, result=[])))

        # 4. Plasmid (PlasmidFinder)
        plasmid_db = args.db / "plasmidfinder_db"
        if plasmid_db.exists():
            tasks.append(task_feature_search("Plasmid_Replicons", plasmid_db, genome_db_prefix, res_dir))
        else:
            tasks.append(asyncio.create_task(asyncio.sleep(0, result=[])))
            
        # 5. Virulence
        virulence_db = args.db / "virulencefinder_db"
        if not virulence_db.exists(): virulence_db = args.db / "Pathogenfinder"
        if virulence_db.exists():
             tasks.append(task_feature_search("Virulence_Factors", virulence_db, genome_db_prefix, res_dir))
        else:
             tasks.append(asyncio.create_task(asyncio.sleep(0, result=[])))

        # 6. Mobile Genetic Elements
        mge_db = args.db / "MEFinder"
        if mge_db.exists():
            tasks.append(task_feature_search("Mobile_Genetic_Elements", mge_db, genome_db_prefix, res_dir))
        else:
            tasks.append(asyncio.create_task(asyncio.sleep(0, result=[])))

        # 실행
        results = await asyncio.gather(*tasks)

        data_collection = {
            "genome_id": genome_id,
            "species_rpoB": results[0],
            "mlst": results[1],
            "amr": results[2],
            "plasmid": results[3],
            "virulence": results[4],
            "mge": results[5]
        }

        write_structured_report(res_dir / "Final_Report.txt", data_collection)

    except Exception as e:
        print(f"[ERROR] 런타임 오류 발생: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
    finally:
        if not args.audit_mode and temp_dir.exists():
            shutil.rmtree(temp_dir)
            print("[SYSTEM] 임시 파일 정리 완료.")
        
    print(f"[SYSTEM] 완료 ({time.time() - start_time:.2f}초)")


# ==============================================================================
# 분석 wrapper
# ==============================================================================

async def make_blast_db(fasta_input: Path, db_prefix: Path):
    if all(db_prefix.with_suffix(s).exists() for s in ['.nin', '.nhr', '.nsq']): return
    cmd = ["makeblastdb", "-in", str(fasta_input), "-dbtype", "nucl", "-out", str(db_prefix), "-parse_seqids"]
    await execute_command(cmd)

async def task_rpoB_wrapper(script_path: Path, genome_fa: Path, db_dir: Path, out_dir: Path) -> Dict:
    print(f"[TASK] rpoB 분석 시작")
    
    abs_script_path = script_path.resolve()
    script_working_dir = abs_script_path.parent
    script_filename = abs_script_path.name
    abs_genome_fa = genome_fa.resolve()
    
    final_out_dir = out_dir / "rpoB_Identification"
    final_out_dir.mkdir(exist_ok=True, parents=True)

    env = os.environ.copy()
    env["RPOB_DB_DIR"] = str(db_dir.resolve())
    ref_core_path = db_dir.parent / "ref" / "rpoB_core.fna"
    if ref_core_path.exists():
        env["RPOB_REF_CORE"] = str(ref_core_path.resolve())
    
    cmd = [sys.executable, script_filename, str(abs_genome_fa)]
    
    proc = await asyncio.create_subprocess_exec(
        *cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE,
        cwd=script_working_dir, env=env
    )
    stdout, stderr = await proc.communicate()
    
    if proc.returncode != 0:
        return {"status": "FAILED", "error": stderr.decode('utf-8', errors='ignore')}

    generated_results_dir = script_working_dir / "results"
    top_hit = "Unknown"
    
    if generated_results_dir.exists():
        for f in generated_results_dir.glob("*"):
            shutil.move(str(f), str(final_out_dir / f.name))
        try: generated_results_dir.rmdir() 
        except: pass

        tsv_path = final_out_dir / "diversity_hits.tsv"
        if tsv_path.exists() and tsv_path.stat().st_size > 0:
            try:
                df = pd.read_csv(tsv_path, sep='\t', header=None)
                if not df.empty and df.shape[1] >= 2:
                    raw = str(df.iloc[0, 1])
                    top_hit = raw.split("|")[0] if "|" in raw else raw
            except: pass

    return {"status": "SUCCESS", "species": top_hit}

async def task_mlst(species: str, genome_db: Path, db_root: Path, out_dir: Path, temp_dir: Path) -> Dict:
    print(f"[TASK] MLST 분석 시작 ({species})")
    species_dir = db_root / "MLST_DB" / species
    
    profile_path = species_dir / f"{species}.txt"
    if not profile_path.exists():
        found = list(species_dir.glob("*.txt"))
        if not found: return {"st": "Unknown (Profile Missing)", "profile": {}}
        profile_path = found[0]

    # 헤더에서 실제 tfa 파일이 있는 Locus만 필터링 (메타데이터 컬럼 제외)
    valid_loci = []
    with open(profile_path, 'r') as f:
        raw_headers = f.readline().strip().split('\t')
        # 첫 컬럼(보통 ST) 제외하고 나머지 중 tfa 파일 있는것만
        potential_loci = raw_headers[1:]
        for l in potential_loci:
            if (species_dir / f"{l}.tfa").exists():
                valid_loci.append(l)
    
    if not valid_loci:
         return {"st": "Unknown (No Valid Loci Found)", "profile": {}}

    print(f"       -> Valid Loci: {valid_loci}") # 디버그용 출력

    probes_fa = temp_dir / "mlst_probes.fasta"
    with open(probes_fa, 'w') as f_out:
        for locus in valid_loci:
            tfa = species_dir / f"{locus}.tfa"
            f_out.write(tfa.read_text())
    
    # 1차 매핑 (Genome vs Probes)
    map_tsv = temp_dir / "mlst_map.tsv"
    await execute_command(["blastn", "-query", str(probes_fa), "-db", str(genome_db), "-out", str(map_tsv), "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", "-perc_identity", "95", "-num_threads", NUM_THREADS])

    # 추출
    extracted_fa = temp_dir / "mlst_extracted.fasta"
    try:
        df = pd.read_csv(map_tsv, sep='\t', names=['q','s','id','len','mis','gap','qs','qe','ss','se','e','bit'])
        best_hits = df.loc[df.groupby('q')['bit'].idxmax()]
    except: return {"st": "Unknown (Detection Fail)", "profile": {}}

    processed = set()
    with open(extracted_fa, 'w') as f_ext:
        for _, row in best_hits.iterrows():
            locus = row['q'].split('_')[0]
            if locus in processed: continue
            s, e = sorted((int(row['ss']), int(row['se'])))
            strand = "plus" if int(row['ss']) < int(row['se']) else "minus"
            succ, seq, _ = await execute_command(["blastdbcmd", "-db", str(genome_db), "-entry", row['s'], "-range", f"{s}-{e}", "-strand", strand])
            if succ:
                f_ext.write(f">{locus}\n{''.join(seq.strip().splitlines()[1:])}\n")
                processed.add(locus)

    # 2차 동정 (Extracted vs All Alleles)
    all_alleles_fa = temp_dir / "mlst_all_alleles.fasta"
    with open(all_alleles_fa, 'w') as f_out:
        for locus in valid_loci:
            f_out.write((species_dir / f"{locus}.tfa").read_text())
            
    allele_db_idx = temp_dir / "mlst_allele_db"
    await make_blast_db(all_alleles_fa, allele_db_idx)

    ident_tsv = temp_dir / "mlst_ident.tsv"
    await execute_command(["blastn", "-query", str(extracted_fa), "-db", str(allele_db_idx), "-out", str(ident_tsv), "-outfmt", "6 qseqid sseqid pident length", "-perc_identity", "100", "-num_threads", NUM_THREADS])

    # 프로파일 매핑
    profile_map = {}
    try:
        df_id = pd.read_csv(ident_tsv, sep='\t', names=['q','s','id','len'])
        df_id = df_id.sort_values('id', ascending=False)
        for locus in valid_loci:
            hit = df_id[df_id['q'] == locus]
            if not hit.empty:
                top = hit.iloc[0]
                match = re.search(r'_(\d+)$', top['s'])
                num = match.group(1) if match else "?"
                profile_map[locus] = num if top['id'] == 100.0 else f"~{num}"
            else: 
                profile_map[locus] = "-"
    except: pass

    # ST 찾기
    st = "Unknown"
    try:
        df_prof = pd.read_csv(profile_path, sep='\t').astype(str)
        q_vec = [profile_map.get(l, "-") for l in valid_loci]
        
        for _, row in df_prof.iterrows():
            # [수정] valid_loci만 비교
            db_vec = [str(row[l]) for l in valid_loci]
            if q_vec == db_vec:
                st = row.iloc[0] # 첫 컬럼이 보통 ST
                break
    except Exception as e:
        print(f"[ERROR] MLST Profile Matching Fail: {e}")

    return {"st": st, "profile": profile_map}

async def task_feature_search(name: str, db_path: Path, genome_db: Path, out_dir: Path) -> List[Dict]:
    print(f"[TASK] {name} 탐색 시작")
    task_dir = out_dir / name
    task_dir.mkdir(exist_ok=True)

    refs = list(db_path.rglob("*.f*a"))
    if not refs: return []
    
    combined_fa = task_dir / "combined_ref.fasta"
    with open(combined_fa, 'w') as f:
        for r in refs: f.write(r.read_text())

    res_tsv = task_dir / "results.tsv"
    await execute_command(["blastn", "-query", str(combined_fa), "-db", str(genome_db), "-out", str(res_tsv), "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp", "-perc_identity", str(DEFAULT_IDENTITY), "-qcov_hsp_perc", str(DEFAULT_COVERAGE), "-num_threads", NUM_THREADS])

    hits = []
    try:
        cols = ['q','s','id','len','mis','gap','qs','qe','ss','se','e','bit','cov']
        df = pd.read_csv(res_tsv, sep='\t', names=cols)
        if not df.empty:
            # 중복 제거 (같은 유전자가 여러 contig에 잡힐 경우 bitscore 높은 것 우선)
            df = df.sort_values('bit', ascending=False).drop_duplicates('q')
            for _, row in df.iterrows():
                hits.append({
                    "Gene": row['q'],
                    "Identity": row['id'],
                    "Coverage": row['cov']
                })
    except: pass
    with open(task_dir / f"{name}.json", 'w') as f: json.dump(hits, f, indent=4)
    return hits

# ==============================================================================
# 리포트 작성
# ==============================================================================

async def execute_command(cmd: List[str], cwd: Path = None) -> dir:
    cmd_str = [str(c) for c in cmd]
    try:
        proc = await asyncio.create_subprocess_exec(
            *cmd_str, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE, cwd=cwd 
        )
        out, err = await proc.communicate()
        return (proc.returncode == 0, out.decode('utf-8', errors='ignore'), err.decode('utf-8', errors='ignore'))
    except Exception as e: return (False, "", str(e))

def find_script_smart(base_path: Path, script_name: str) -> Path:
    if (base_path / script_name).exists(): return (base_path / script_name).resolve()
    found = list(base_path.rglob(script_name))
    return found[0].resolve() if found else None

def write_structured_report(path: Path, data: Dict):
    with open(path, 'w') as f:
        f.write("REPORT\n\n")
        
        # 1. Species
        species_val = data['species_rpoB'].get('species', 'Unknown')
        f.write(f"Species: {species_val}\n\n")

        # 2. Molecular epidemiology
        f.write("Molecular epidemiology\n")
        mlst_info = data['mlst']
        st_val = mlst_info.get('st', 'Unknown')
        f.write(f"  MLST: {st_val}\n")
        # [수정] 상세 Allele 정보 출력 (디버깅 및 정보 제공용)
        profile = mlst_info.get('profile', {})
        if profile:
            prof_str = ", ".join([f"{k}({v})" for k,v in profile.items()])
            f.write(f"  Alleles: {prof_str}\n")
        
        # 3. AMR
        f.write("Antimicrobial resistance determinants\n")
        acquired_genes = [f"{h['Gene']} ({h['Identity']}%)" for h in data['amr']] if data['amr'] else []
        
        f.write("  acquired genes\n")
        if acquired_genes:
            for gene in acquired_genes: f.write(f"    {gene}\n")
        else: f.write("    None\n")
        f.write("  SNP\n    None\n")

        # 4. Virulence
        f.write("Virulence factors\n")
        vir = [f"{h['Gene']} ({h['Identity']}%)" for h in data.get('virulence', [])]
        if vir:
            for v in vir: f.write(f"  {v}\n")
        else: f.write("  None\n")

        # 5. MGE
        f.write("Mobile genetic elements\n")
        f.write("  plasmid\n")
        plasmids = [f"{h['Gene']} ({h['Identity']}%)" for h in data.get('plasmid', [])]
        if plasmids:
            for p in plasmids: f.write(f"    {p}\n")
        else: f.write("    None\n")
            
        f.write("  mobile genetic elements\n")
        mges = [f"{h['Gene']} ({h['Identity']}%)" for h in data.get('mge', [])]
        if mges:
            for m in mges: f.write(f"    {m}\n")
        else: f.write("    None\n")

if __name__ == "__main__":
    asyncio.run(main())