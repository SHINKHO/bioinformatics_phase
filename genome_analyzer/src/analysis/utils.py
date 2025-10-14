# src/analysis/utils.py
import subprocess
from pathlib import Path
from Bio import SeqIO
# [수정] MLST_DB_ROOT 대신 DATABASE_ROOT를 import 합니다.
from config import DATABASE_ROOT

def check_dependencies():
    """Checks if required command-line tools are available."""
    dependencies = ["blastn", "makeblastdb", "blastdbcmd"]
    for dep in dependencies:
        if subprocess.run(["which", dep], capture_output=True, text=True).returncode != 0:
            raise RuntimeError(
                f"Dependency '{dep}' not found in PATH. "
                "Please install NCBI BLAST+ and ensure it's in your system's PATH."
            )

def setup_mlst_parameters(genome_file: Path, logger):
    """
    Identifies species from the genome file and sets up MLST-related paths and data.
    """
    # 1. Identify species from the first record's description
    try:
        first_record = next(SeqIO.parse(genome_file, "fasta"))
        header = first_record.description.lower()
        logger.log_step("MLST", "1_1_Read_Genome_Header", f"Reading genome file header: {header}")
    except StopIteration:
        raise ValueError(f"FASTA file '{genome_file}' is empty or not in a valid format.")

    species_match = None
    potential_species = ["klebsiella", "escherichia", "salmonella", "staphylococcus", "streptococcus", "enterococcus"]
    for sp in potential_species:
        if sp in header:
            species_match = sp
            logger.log_step("MLST", "1_2_Species_Name_Extraction", f"Extracted species name: {species_match}")
            break

    if not species_match:
        raise ValueError(
            "Could not automatically determine the species from the FASTA header. "
            "Please ensure the header contains a recognizable species name (e.g., 'Klebsiella pneumoniae')."
        )

    # 2. [수정] DATABASE_ROOT를 기준으로 동적으로 경로를 생성합니다.
    species_db_dir = DATABASE_ROOT / "MLST_DB" / species_match
    if not species_db_dir.is_dir():
        raise FileNotFoundError(f"No MLST database found for species '{species_match}' at '{species_db_dir}'")

    profile_file = next(species_db_dir.glob("*.txt"), None)
    if not profile_file:
        raise FileNotFoundError(f"MLST profile file (.txt) not found in '{species_db_dir}'")

    # 3. Determine the order of loci from the profile file header
    with open(profile_file, 'r') as f:
        header_line = f.readline().strip()
        loci_order = header_line.split('\t')[1:] # Skip the 'ST' column

    return {
        "species": species_match,
        "gene_dir": species_db_dir,
        "profile_file": profile_file,
        "loci_order": loci_order
    }