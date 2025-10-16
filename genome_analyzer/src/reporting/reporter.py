"""
Final Report Generator

This module is responsible for parsing all analysis results and generating
the final, human-readable 'One Page ME Report'. It consolidates data from
MLST, AMR, plasmid, and MGE analyses into a single text file.
"""
from pathlib import Path
import pandas as pd

# Project-level module imports
from config import ANALYSES_TO_RUN

def create_final_report(results_data: dict, results_dir: Path, genome_name: str):
    """
    Generates the final summary report from the collected analysis data.

    This function orchestrates the creation of the final report file. It writes
    a header with sample information, followed by dedicated sections for MLST
    and each of the other analyses defined in the configuration.

    Args:
        results_data (dict): The main dictionary containing all analysis results.
                             It holds MLST data under the 'mlst' key.
        results_dir (Path): The top-level directory where all results, including this
                            report, will be saved.
        genome_name (str): The identifier for the genome, extracted from the FASTA file.
    """
    # Step 1: Define the path for the final report file.
    report_path = results_dir / "Final_ME_Report.txt"
    
    # Step 2: Extract MLST results and parameters from the main results dictionary.
    mlst_params = results_data.get('mlst_params', {})
    mlst_results = results_data.get('mlst', {})
    
    # Step 3: Open the report file for writing.
    with open(report_path, "w") as f:
        # --- Section A: Header ---
        # Step 4: Write the main title and basic sample information.
        f.write("========== One Page ME Report ==========\n\n")
        f.write(f"■ ID: {genome_name}\n")
        f.write(f"■ Species: {mlst_params.get('species', 'Unknown')}\n\n")

        # --- Section B: MLST Results ---
        # Step 5: Write the Molecular Epidemiology section, focusing on MLST.
        f.write("--- Molecular Epidemiology ---\n")
        st = mlst_results.get('st', 'Not determined')
        profile = mlst_results.get('alleles', {})
        profile_str = ", ".join([f"{locus}-{num}" for locus, num in profile.items()])
        f.write(f"▶ MLST: {st}\n")
        f.write(f"  Allele Profile: {profile_str}\n\n")

        # --- Section C: Other BLAST Analyses ---
        # Step 6: Loop through all analyses defined in the config to report their results.
        col_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        
        # A dictionary to group analyses under common headers.
        analysis_groups = {
            "Antimicrobial Resistance Determinants": ["Antimicrobial_Resistance"],
            "Mobile Genetic Elements": ["Plasmid_Replicons", "Mobile_Genetic_Elements"]
        }

        # Step 7: Process each analysis group.
        for group_title, analysis_names in analysis_groups.items():
            f.write(f"--- {group_title} ---\n")
            found_any_in_group = False

            # Step 8: Process each analysis within the group.
            for analysis_name in analysis_names:
                # Find the corresponding db_folder from ANALYSES_TO_RUN
                db_folder = next((db for db, name in ANALYSES_TO_RUN.items() if name == analysis_name), None)
                if not db_folder:
                    continue

                # Determine the subtitle for the analysis.
                if "Resistance" in analysis_name:
                    subtitle = "Acquired Genes"
                elif "Replicons" in analysis_name:
                    subtitle = "Plasmid Replicons"
                else:
                    subtitle = "Other MGEs"
                f.write(f"▶ {subtitle}:\n")

                # Step 9: Read the BLAST result file for the current analysis.
                result_file = results_dir / analysis_name / "blast_results.tsv"
                try:
                    df = pd.read_csv(result_file, sep='\t', names=col_names)
                    if df.empty:
                        # Step 10a: If no hits, report that.
                        f.write("  - No significant hits found.\n")
                    else:
                        # Step 10b: If hits are found, get the best hit for each query sequence.
                        best_hits = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]
                        for _, row in best_hits.iterrows():
                            f.write(f"  - {row['qseqid']:<25} (Identity: {row['pident']:.2f}%, Contig: {row['sseqid']})\n")
                        found_any_in_group = True
                except (FileNotFoundError, pd.errors.EmptyDataError):
                    # Step 10c: Handle cases where the result file doesn't exist or is empty.
                    f.write("  - No significant hits found.\n")
            f.write("\n")

    # Step 11: Print a final confirmation message to the console.
    print(f"✅ Final report created at: {report_path}")