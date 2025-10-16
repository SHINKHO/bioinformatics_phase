"""
Utility Functions for the Analysis Pipeline

This module provides helper functions that support the main analysis workflow.
This includes checking for external dependencies and preparing parameters for
specific analyses like MLST.
"""
import subprocess
from pathlib import Path
from Bio import SeqIO
from typing import Dict, Any, List

# Project-level module imports
from config import DATABASE_ROOT
from logger import Logger


def check_dependencies():
    """
    Checks if required command-line tools are installed and in the system's PATH.

    This function iterates through a list of essential NCBI BLAST+ tools and uses
    the `which` command to verify their existence. If any tool is not found, it
    indicates a setup problem.

    Raises:
        RuntimeError: If any of the dependency tools are not found in the PATH.
    """
    # Step 1: Define the list of required command-line dependencies.
    dependencies = ["blastn", "makeblastdb", "blastdbcmd", "prodigal", "python", "diamond"]
    
    # Step 2: Loop through each dependency.
    for dep in dependencies:
        # Step 3: Use `which` to check if the command exists in the system's PATH.
        # A non-zero return code means the command was not found.
        if subprocess.run(["which", dep], capture_output=True, text=True).returncode != 0:
            # Step 4: If not found, raise an error with an informative message.
            raise RuntimeError(
                f"Dependency '{dep}' not found in PATH. "
                "Please install NCBI BLAST+ and ensure it's in your system's PATH."
            )


def setup_mlst_parameters(genome_file: Path, logger: Logger) -> Dict[str, Any]:
    """
    Identifies the species from the genome folder structure and prepares MLST-specific parameters.

    This function assumes the genome file is located in a folder structure like:
    {genome_id}/{species}/genome_file.fasta
    It extracts the species name from the parent folder and uses it to locate the
    corresponding MLST database, profile file, and determine the correct order of loci.

    Args:
        genome_file (Path): The path to the input genome FASTA file.
        logger (Logger): The logger instance for detailed step-logging.

    Returns:
        dict: A dictionary containing all necessary parameters for the MLST workflow:
              - "species" (str): The identified species name.
              - "gene_dir" (Path): The path to the species-specific MLST database directory.
              - "profile_file" (Path): The path to the MLST profile definition file.
              - "loci_order" (list): A list of locus names in the correct order.
              - "genome_id" (str): The genome identifier from the folder structure.

    Raises:
        ValueError: If the folder structure is invalid or the species cannot be determined.
        FileNotFoundError: If the MLST database or profile file for the species is not found.
    """
    # Step 1: Extract species and genome_id from folder structure.
    # Expected structure: {genome_id}/{species}/genome_file.fasta
    try:
        parent_dir = genome_file.parent
        species_dir = parent_dir.name
        genome_id = parent_dir.parent.name
        
        logger.log_step("MLST", "1_1_Read_Folder_Structure",
                       f"Extracted from folder structure: genome_id='{genome_id}', species='{species_dir}'")
    except Exception as e:
        raise ValueError(
            f"Invalid genome folder structure. Expected format: {{genome_id}}/{{species}}/genome_file.fasta. "
            f"Got: '{genome_file}'. Error: {e}"
        )

    # Step 2: Validate that the species exists in the MLST database.
    mlst_db_path = DATABASE_ROOT / "MLST_DB"
    potential_species = [d.name for d in mlst_db_path.iterdir() if d.is_dir()]
    
    if species_dir not in potential_species:
        raise ValueError(
            f"Species '{species_dir}' not found in MLST database. "
            f"Available species: {', '.join(potential_species)}"
        )

    # Step 3: Construct the path to the species-specific MLST database.
    species_db_dir = mlst_db_path / species_dir
    logger.log_step("MLST", "1_2_Species_Database_Found", f"MLST database found for species '{species_dir}' at '{species_db_dir}'")

    # Step 4: Find the MLST profile file (e.g., klebsiella.txt) in the database directory.
    profile_file = next(species_db_dir.glob("*.txt"), None)
    if not profile_file:
        raise FileNotFoundError(f"MLST profile file (.txt) not found in '{species_db_dir}'")

    # Step 5: Determine the correct order of loci from the header of the profile file.
    # Only read the first few lines to avoid reading large files
    with open(profile_file, 'r') as f:
        header_line = f.readline().strip()
        # The loci are the columns after the first 'ST' column.
        loci_order = header_line.split('\t')[1:]

    # Step 6: Return all the prepared parameters in a dictionary.
    return {
        "species": species_dir,
        "gene_dir": species_db_dir,
        "profile_file": profile_file,
        "loci_order": loci_order,
        "genome_id": genome_id
    }

def get_genome_name(genome_file: Path) -> str:
    """
    Extracts a descriptive name from the genome FASTA file.

    This function reads the first record of a FASTA file and returns its ID,
    which is typically a unique accession number or identifier.

    Args:
        genome_file (Path): The path to the input genome FASTA file.

    Returns:
        str: The identifier of the first sequence record.

    Raises:
        ValueError: If the FASTA file is empty or not in a valid format.
    """
    try:
        first_record = next(SeqIO.parse(genome_file, "fasta"))
        return first_record.id
    except StopIteration:
        raise ValueError(f"FASTA file '{genome_file}' is empty or not in a valid format.")