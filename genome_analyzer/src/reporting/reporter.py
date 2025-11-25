# 최종 분석 결과를 파싱하고 사람이 읽을 수 있는 보고서를 생성.
from pathlib import Path
import pandas as pd

# Project-level module imports
from config import ANALYSES_TO_RUN

def create_final_report(results_data: dict, results_dir: Path, genome_name: str):
    # 수집된 분석 데이터로부터 최종 요약 보고서를 생성. In: results_data(dict), results_dir(Path), genome_name(str) / Out: None
    report_path = results_dir / "Final_ME_Report.txt"
    
    mlst_params = results_data.get('mlst_params', {})
    mlst_results = results_data.get('mlst', {})
    
    with open(report_path, "w") as f:
        # 헤더
        f.write("========== One Page ME Report ==========\n\n")
        f.write(f"■ ID: {genome_name}\n")
        f.write(f"■ Species: {mlst_params.get('species', 'Unknown')}\n\n")

        # MLST 결과
        f.write("--- Molecular Epidemiology ---\n")
        st = mlst_results.get('st', 'Not determined')
        profile = mlst_results.get('alleles', {})
        profile_str = ", ".join([f"{locus}-{num}" for locus, num in profile.items()])
        f.write(f"▶ MLST: {st}\n")
        f.write(f"  Allele Profile: {profile_str}\n\n")

        # 기타 BLAST 분석 결과
        col_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        
        analysis_groups = {
            "Antimicrobial Resistance Determinants": ["Antimicrobial_Resistance"],
            "Mobile Genetic Elements": ["Plasmid_Replicons", "Mobile_Genetic_Elements"]
        }

        for group_title, analysis_names in analysis_groups.items():
            f.write(f"--- {group_title} ---\n")
            found_any_in_group = False

            for analysis_name in analysis_names:
                db_folder = next((db for db, name in ANALYSES_TO_RUN.items() if name == analysis_name), None)
                if not db_folder:
                    continue

                if "Resistance" in analysis_name:
                    subtitle = "Acquired Genes"
                elif "Replicons" in analysis_name:
                    subtitle = "Plasmid Replicons"
                else:
                    subtitle = "Other MGEs"
                f.write(f"▶ {subtitle}:\n")

                result_file = results_dir / analysis_name / "blast_results.tsv"
                try:
                    df = pd.read_csv(result_file, sep='\t', names=col_names)
                    if df.empty:
                        f.write("  - No significant hits found.\n")
                    else:
                        best_hits = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]
                        for _, row in best_hits.iterrows():
                            f.write(f"  - {row['qseqid']:<25} (Identity: {row['pident']:.2f}%, Contig: {row['sseqid']})\n")
                        found_any_in_group = True
                except (FileNotFoundError, pd.errors.EmptyDataError):
                    f.write("  - No significant hits found.\n")
            f.write("\n")

    print(f"✅ Final report created at: {report_path}")