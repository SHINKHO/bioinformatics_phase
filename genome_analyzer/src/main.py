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
import json
from pathlib import Path

from config import DEFAULT_RESULTS_DIR, DATABASE_ROOT
from pipeline_service import PipelineService

def main():
    """
    Parses command-line arguments and starts the analysis pipeline.
    
    This function defines the user-facing command-line interface (CLI) and
    passes the parsed arguments to the PipelineService to begin the workflow.
    """
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
    
    parser.add_argument(
        "-p", "--pipeline",
        type=Path,
        default=Path("pipelines/default.json"),
        help="Path to the pipeline configuration JSON file.\n(default: pipelines/default.json)"
    )
    
    args = parser.parse_args()

    # Load the pipeline configuration from JSON
    with open(args.pipeline, "r") as f:
        pipeline_config = json.load(f)

    # Prepare the initial context for the service
    initial_context = {
        "initial_genome_file": args.genome_file,
        "base_results_dir": args.output,
        "verbose": args.verbose,
        "database_root": DATABASE_ROOT
    }
    
    # Instantiate and run the pipeline service
    service = PipelineService(pipeline_config, initial_context)
    asyncio.run(service.run_analysis())

if __name__ == "__main__":
    # This block ensures that main() is called only when the script is executed directly.
    main()