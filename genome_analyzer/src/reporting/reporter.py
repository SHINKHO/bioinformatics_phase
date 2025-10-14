"""
Report Generator

This module is responsible for parsing all analysis results and generating
the final, human-readable 'One Page ME Report'.
"""
from pathlib import Path
import pandas as pd
from config import ANALYSES_TO_RUN

def create_final_report(results_data: dict, results_dir: Path, genome_file: Path):
    """
    Generates the final summary report from the collected analysis data.
    """
    report_path = results_dir / "Final_ME_Report.txt"
    mlst_params = results_data.get('mlst_params', {})
    mlst_results = results_data.get('mlst', {})
    
    with open(report_path, "w") as f:
        f.write("========== One Page ME Report ==========\n\n")
        f.write(f"■ ID: {genome_file.name}\n")
        f.write(f"■ Species: {mlst_params.get('species', 'Unknown')}\n\n")

        # --- MLST results ---
        f.write("--- Molecular Epidemiology ---\n")
        st = mlst_results.get('ST', 'Not determined')
        profile = mlst_results.get('Profile', {})
        profile_str = ", ".join([f"{locus}-{num}" for locus, num in profile.items()])
        f.write(f"▶ MLST: {st}\n")
        f.write(f"  Allele Profile: {profile_str}\n\n")

        # --- Other analysis results ---
        col_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        
        for db_folder, analysis_name in ANALYSES_TO_RUN.items():
            title = analysis_name.replace("_", " ")
            header_written = False

            if "Resistance" in title:
                f.write("--- Antimicrobial Resistance Determinants ---\n▶ Acquired Genes:\n")
                header_written = True
            elif "Plasmid" in title:
                f.write("--- Mobile Genetic Elements ---\n▶ Plasmid Replicons:\n")
                header_written = True
            elif "Mobile" in title:
                f.write("▶ Other MGEs:\n")
                header_written = True

            result_file = results_dir / analysis_name / "blast_results.tsv"
            try:
                df = pd.read_csv(result_file, sep='\t', names=col_names)
                if df.empty:
                    f.write("  - No significant hits found.\n")
                else:
                    best_hits = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]
                    for _, row in best_hits.iterrows():
                        f.write(f"  - {row['qseqid']:<25} (Identity: {row['pident']:.2f}%, Contig: {row['sseqid']})\n")
            except (FileNotFoundError, pd.errors.EmptyDataError):
                f.write("  - No significant hits found.\n")
            
            if header_written:
                f.write("\n")

    print(f"✅ Final report created at: {report_path}")
