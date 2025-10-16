"""
Genome Analysis Pipeline - Main Runner

This script serves as the main entry point for the genome analysis pipeline.
Its primary responsibilities are:
1. Parsing command-line arguments.
2. Instantiating the AnalysisManager.
3. Initiating the asynchronous analysis workflow.
"""
import argparse
import asyncio
from pathlib import Path

# Project-level module imports
from config import DEFAULT_RESULTS_DIR
from analysis.manager import AnalysisManager

def main():
    """
    Parses command-line arguments and starts the analysis pipeline.
    
    This function defines the user-facing command-line interface (CLI) and
    passes the parsed arguments to the AnalysisManager to begin the workflow.
    """
    # Step 1: Create the argument parser with a description.
    parser = argparse.ArgumentParser(
        description="A comprehensive genome analysis pipeline for molecular epidemiology.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Step 2: Define the command-line arguments.
    # The input genome file (required positional argument).
    parser.add_argument(
        "genome_file",
        type=Path,
        help="Path to the input genome file in FASTA format."
    )
    # The output directory (optional, with a default).
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help=f"Path to the output directory.\n(default: {DEFAULT_RESULTS_DIR})"
    )
    # The verbose flag (optional).
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose mode to see detailed progress logs."
    )
    
    # Step 3: Parse the arguments provided by the user.
    args = parser.parse_args()
    
    # Step 4: Instantiate the main orchestrator, AnalysisManager, with the parsed arguments.
    manager = AnalysisManager(
        genome_file=args.genome_file,
        results_dir=args.output,
        verbose=args.verbose
    )
    
    # Step 5: Run the main asynchronous pipeline.
    # asyncio.run() starts the event loop and runs the coroutine until it completes.
    asyncio.run(manager.run_pipeline())

if __name__ == "__main__":
    # This block ensures that main() is called only when the script is executed directly.
    main()