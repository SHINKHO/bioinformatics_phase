"""
Main Analysis Orchestrator

This module contains the AnalysisManager class, which is the main orchestrator for
the entire genome analysis pipeline.
"""
import asyncio
import time
from pathlib import Path
import shutil

# Local (sibling) module imports
from . import utils
from . import blast_runner
from .handler import (AnalysisContext, MLSTHandler, PathogenFinder2Handler, 
                          StandardAnalysisHandler, AMRHandler)

# Project-level module imports
from reporting import reporter
from config import ANALYSES_TO_RUN, BLAST_DB_DIR
from logger import Logger

class AnalysisManager:
    """
    Orchestrates the genome analysis workflow from start to finish.
    
    This class is responsible for:
    - Setting up the environment.
    - Creating a BLAST database from the input genome.
    - Dispatching analysis tasks to a chain of handlers.
    - Generating the final report.
    - Cleaning up temporary files.
    """
    def __init__(self, genome_file: Path, results_dir: Path, verbose: bool = False):
        """
        Initializes the AnalysisManager.
        
        Args:
            genome_file (Path): Path to the input genome file in FASTA format.
            results_dir (Path): Path to the output directory for results.
            verbose (bool): Flag to enable verbose console logging.
        """
        self.genome_file = genome_file
        self.base_results_dir = results_dir # Store base results dir
        self.results_dir = results_dir # To be updated
        self.temp_dir = results_dir / "temp" # To be updated
        self.verbose = verbose
        self.results_data = {}
        
        # Set up logging
        self.base_logs_dir = Path.cwd() / "logs"
        self.logs_dir = self.base_logs_dir
        self.logger = Logger(self.logs_dir)

    def _log(self, message: str, level: str = "INFO"):
        """
        Prints a log message to the console if in verbose mode.
        
        Args:
            message (str): The message to print.
            level (str): The log level (e.g., "INFO", "WARN").
        """
        if self.verbose:
            timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
            print(f"[{timestamp} - {level}] {message}")

    async def run_pipeline(self):
        """
        Executes the entire analysis pipeline from start to finish.
        
        This is the main entry point method. It orchestrates all steps:
        setup, DB creation, concurrent analysis, reporting, and cleanup.
        """
        start_time = time.time()
        print("===== Genome Analysis Pipeline Start =====")
        try:
            # --- Step 1: Pre-flight checks and setup ---
            self._log("Step 1: Pre-flight checks and setup.")
            self.logger.log_step("Pipeline", "1_Pre-flight_Checks", "Starting pre-flight checks and setup.")
            utils.check_dependencies()

            # --- Step 2: Identify species and set up paths ---
            self._log("Step 2: Identifying species and setting up paths.")
            mlst_params = utils.setup_mlst_parameters(self.genome_file, self.logger)
            self.results_data['mlst_params'] = mlst_params
            genome_id = mlst_params['genome_id']
            species = mlst_params['species']
            self._log(f"Species '{species}' identified for MLST from folder structure.")
            self.logger.log_step("Pipeline", "2_Species_Identification", f"Species '{species}' identified from folder structure.")

            # Define output directories based on genome ID and species
            self.results_dir = self.base_results_dir / genome_id / species
            self.temp_dir = self.results_dir / "temp"
            self.logs_dir = self.base_logs_dir / genome_id / species
            blast_db_dir = BLAST_DB_DIR / genome_id / species

            # Create directories
            self.results_dir.mkdir(parents=True, exist_ok=True)
            self.temp_dir.mkdir(exist_ok=True)
            blast_db_dir.mkdir(parents=True, exist_ok=True)
            
            # Re-initialize logger with the new path
            self.logger = Logger(self.logs_dir)

            # --- Step 3: Create BLAST database for the input genome ---
            self._log("Step 3: Creating BLAST database for the input genome.")
            self.logger.log_step("Pipeline", "3_Create_BLAST_DB", "Creating BLAST database for the input genome.")
            genome_db_path = await blast_runner.create_blast_db_async(self.genome_file, blast_db_dir)
            self._log(f"Genome BLAST DB created at '{genome_db_path}'.")
            self.logger.log_step("Pipeline", "4_BLAST_DB_Created", f"Genome BLAST DB created at '{genome_db_path}'.")

            # --- Step 4: Run all analysis tasks concurrently ---
            self._log("Step 4: Running all analysis tasks concurrently.")
            self.logger.log_step("Pipeline", "5_Run_Concurrent_Analyses", "Running all analysis tasks concurrently.")
            
            # 4a. Prepare context for handlers
            context = AnalysisContext(
                genome_db_path=genome_db_path,
                results_dir=self.results_dir,
                temp_dir=self.temp_dir,
                logger=self.logger,
                verbose=self.verbose,
                results_data=self.results_data,
                genome_id=genome_id,
                species=species
            )

            # 4b. Build the chain of responsibility
            # The chain is: MLSTHandler -> AMRHandler -> PathogenFinder2Handler -> StandardAnalysisHandler
            standard_handler = StandardAnalysisHandler(context)
            pathogen_handler = PathogenFinder2Handler(context)
            amr_handler = AMRHandler(context)
            analysis_chain = MLSTHandler(context)
            analysis_chain.set_next(amr_handler).set_next(pathogen_handler).set_next(standard_handler)

            # 4c. Dispatch all analyses to the handler chain
            tasks = []
            for db_folder, analysis_name in ANALYSES_TO_RUN.items():
                # Pass analysis-specific params to the appropriate analysis
                if analysis_name == "MLST":
                    params = mlst_params
                elif analysis_name == "Pathogen_Finder2":
                    # PathogenFinder2-specific parameters
                    params = {
                        "database_dir": str(Path.cwd() / "database" / "Pathogenfinder"),
                        "output_dir": str(self.results_dir / "Pathogen_Finder2"),
                        "genome_path": str(genome_db_path)
                    }
                else:
                    params = {}
                
                # The handle method returns a ready-to-run asyncio.Task
                task = await analysis_chain.handle(
                    analysis_name=analysis_name,
                    db_folder=db_folder,
                    params=params
                )
                if task:
                    tasks.append(task)
            
            # 4d. Run all created tasks concurrently
            await asyncio.gather(*tasks)
            self._log("All analysis tasks completed.")
            self.logger.log_step("Pipeline", "6_Concurrent_Analyses_Complete", "All analysis tasks completed.")

            # --- Step 5: Generate final report ---
            self._log("Step 5: Generating final report.")
            self.logger.log_step("Pipeline", "7_Generate_Report", "Generating final report.")
            genome_name = utils.get_genome_name(self.genome_file)
            reporter.create_final_report(self.results_data, self.results_dir, genome_name)
            self._log("Final report generated.")
            self.logger.log_step("Pipeline", "8_Report_Generated", "Final report generated.")

        except (ValueError, FileNotFoundError, RuntimeError, Exception) as e:
            # --- Error Handling ---
            print(f"\n❌ PIPELINE FAILED: An error occurred.\n  -> {e}")
            self.logger.log_step("Pipeline", "9_Pipeline_Failed", f"PIPELINE FAILED: An error occurred.\n  -> {e}")
        
        finally:
            # --- Step 6: Cleanup ---
            if self.temp_dir.exists():
                shutil.rmtree(self.temp_dir)
        
        end_time = time.time()
        print(f"\n🎉 Analysis complete in {end_time - start_time:.2f} seconds.")
        print(f"   Results are in '{self.results_dir}'")
        print("==========================================")
