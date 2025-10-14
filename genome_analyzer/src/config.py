# src/config.py
from pathlib import Path

# --- Project Root ---
PROJECT_ROOT = Path.cwd()

# --- Main Directories ---
DATABASE_ROOT = PROJECT_ROOT / "database"
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "analysis_results"

# --- BLAST Database ---
BLAST_DB_DIR = PROJECT_ROOT / "blast_db_output"
BLAST_DB_DIR.mkdir(exist_ok=True)

# --- Analyses to Run ---
# Dictionary mapping database folders to the desired output analysis names.
# MLST analysis is now included here for consistent processing.
ANALYSES_TO_RUN = {
    "MLST_DB": "MLST", # <-- MLST 분석 추가
    "resfinder_db": "Antimicrobial_Resistance",
    "plasmidfinder_db": "Plasmid_Replicons",
    "mefinder_db": "Mobile_Genetic_Elements",
}