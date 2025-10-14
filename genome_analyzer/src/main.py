"""
Genome Analysis Pipeline - Main Runner

This is the main entry point for the analysis pipeline.
It parses command-line arguments and initiates the analysis
process using the AnalysisManager.
"""
import argparse
import asyncio
from pathlib import Path
from config import DEFAULT_RESULTS_DIR
from analysis.manager import AnalysisManager

def main():
    """Parses arguments and starts the pipeline."""
    parser = argparse.ArgumentParser(
        description="A comprehensive genome analysis pipeline for molecular epidemiology.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "genome_file",
        type=Path,
        help="Path to the input genome file in FASTA format."
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help=f"Path to the output directory.\n(default: {DEFAULT_RESULTS_DIR})"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose mode to see detailed progress logs."
    )
    
    args = parser.parse_args()
    
    # Instantiate and run the manager
    manager = AnalysisManager(
        genome_file=args.genome_file,
        results_dir=args.output,
        verbose=args.verbose
    )
    
    # Run the asynchronous pipeline
    asyncio.run(manager.run_pipeline())

if __name__ == "__main__":
    main()
