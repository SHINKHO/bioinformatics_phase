"""
Core Pipeline Service

This module provides the main PipelineService class that orchestrates the
data-driven analysis workflow.
"""
import asyncio
import importlib
from typing import Dict, Any
from pathlib import Path

from analysis.context import WorkflowContext
from analysis import utils, blast_runner
from logger import Logger

class PipelineService:
    """
    The core, decoupled pipeline execution engine.
    
    This service reads a pipeline configuration, manages the workflow context,
    and executes each analysis step in sequence. It is designed to be
    framework-agnostic and can be invoked from a CLI, web app, or any other
    trigger.
    """
    def __init__(self, pipeline_config: Dict[str, Any], initial_context: Dict[str, Any]):
        """
        Initializes the PipelineService.
        
        Args:
            pipeline_config (Dict[str, Any]): The pipeline definition, typically loaded from a JSON file.
            initial_context (Dict[str, Any]): A dictionary of initial values for the workflow context
                                             (e.g., input files, root directories).
        """
        self.pipeline_config = pipeline_config
        self.context = WorkflowContext(initial_context)
        self.results = {}

    async def run_analysis(self):
        """
        Runs the entire analysis pipeline as defined by the configuration.
        """
        print(f"PipelineService: Starting pipeline '{self.pipeline_config.get('pipeline_name', 'Untitled')}'...")
        
        # --- Step 1: Pre-flight checks and setup (from old AnalysisManager) ---
        utils.check_dependencies()
        
        initial_genome_file = Path(self.context.get("initial_genome_file"))
        base_results_dir = Path(self.context.get("base_results_dir"))
        base_logs_dir = Path.cwd() / "logs"

        # Identify species and get initial params
        # Note: This assumes the MLST folder structure for species identification.
        # This could be made more generic in a future refactoring.
        mlst_params = utils.setup_mlst_parameters(initial_genome_file, Logger(base_logs_dir))
        
        # Populate context with setup info
        self.context["genome_id"] = mlst_params['genome_id']
        self.context["species"] = mlst_params['species']
        self.context["results_dir"] = base_results_dir / self.context["genome_id"] / self.context["species"]
        self.context["temp_dir"] = self.context["results_dir"] / "temp"
        self.context["logs_dir"] = base_logs_dir / self.context["genome_id"] / self.context["species"]
        
        # Create directories
        self.context["results_dir"].mkdir(parents=True, exist_ok=True)
        self.context["temp_dir"].mkdir(exist_ok=True)
        
        # Re-initialize logger with the specific run path
        logger = Logger(self.context["logs_dir"])
        self.context["logger"] = logger
        
        # Create BLAST database for the input genome
        blast_db_output_dir = Path.cwd() / "blast_db_output" / self.context["genome_id"] / self.context["species"]
        blast_db_output_dir.mkdir(parents=True, exist_ok=True)
        genome_db_path = await blast_runner.create_blast_db_async(initial_genome_file, blast_db_output_dir)
        self.context["genome_db_path"] = genome_db_path
        
        # Initialize results data dict in context
        self.context["results_data"] = {}

        # --- Step 2: Execute pipeline steps concurrently ---
        tasks = []
        steps = self.pipeline_config.get("steps", [])
        for i, step_config in enumerate(steps):
            step_name = step_config.get("name", f"Unnamed Step {i+1}")
            
            if not step_config.get("enabled", True):
                logger.log_step("Pipeline", f"Skip_{step_name}", f"Skipping disabled step: {step_name}")
                continue

            logger.log_step("Pipeline", f"Start_{step_name}", f"Initiating step: {step_name}")
            
            module_name = step_config.get("module")
            if not module_name:
                logger.log_step("Pipeline", f"Error_{step_name}", f"ERROR: 'module' not defined for step: {step_name}")
                continue

            try:
                # Dynamically import the handler class
                module_path, class_name = module_name.rsplit('.', 1)
                handler_module = importlib.import_module(module_path)
                handler_class = getattr(handler_module, class_name)
                
                # Resolve I/O paths
                resolved_inputs = {k: self.context.resolve_path(v) for k, v in step_config.get("inputs", {}).items()}
                resolved_outputs = {k: self.context.resolve_path(v) for k, v in step_config.get("outputs", {}).items()}
                parameters = step_config.get("parameters", {})

                # Instantiate handler and create execution task
                handler_instance = handler_class(self.context, step_config)
                task = handler_instance.execute()
                tasks.append(task)

            except (ImportError, AttributeError, ValueError) as e:
                logger.log_step("Pipeline", f"Error_{step_name}", f"ERROR: Could not load or run step {step_name}: {e}")
            
        # Run all created tasks concurrently
        await asyncio.gather(*tasks)
        
        # --- Step 3: Final reporting ---
        from reporting import reporter
        from analysis import utils as analysis_utils

        genome_name = analysis_utils.get_genome_name(initial_genome_file)
        steps = self.pipeline_config.get("steps", [])
        self.context["results_data"]["species"] = self.context.get("species")
        reporter.create_final_report(self.context.get("results_data"), self.context.get("results_dir"), genome_name, steps)
        
        logger.log_step("Pipeline", "Complete", "Pipeline execution complete.")
        print("PipelineService: Analysis complete.")
        return self.context.get("results_data")

