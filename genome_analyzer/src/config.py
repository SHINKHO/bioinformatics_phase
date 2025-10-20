"""
Central Configuration for the Pipeline

This module centralizes all the main configuration variables for the project.
This includes file paths, directory names, and the crucial dictionary that
defines which analyses to run.
"""
from pathlib import Path

# --- Project Root ---
# Defines the root directory of the project. `Path.cwd()` assumes the
# script is run from the project's root directory.
PROJECT_ROOT = Path.cwd()

# --- Main Directories ---
# Defines the primary data and result directories, relative to the project root.
DATABASE_ROOT = PROJECT_ROOT / "database"
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "analysis_results"

# --- BLAST Database Directory ---
# This directory is used to store the pre-formatted BLAST databases created
# from the input genomes. This avoids re-creating the database on every run
# for the same genome.
BLAST_DB_DIR = PROJECT_ROOT / "blast_db_output"
BLAST_DB_DIR.mkdir(exist_ok=True)


