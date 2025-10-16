
import asyncio
import subprocess
from pathlib import Path
import pandas as pd
import json

from .base import AnalysisHandler
from analysis import pathogen_runner

class PathogenFinder2Handler(AnalysisHandler):
    """
    A concrete handler for the PathogenFinder2 workflow.
    
    This handler checks if the requested analysis is "Pathogen_Finder2". If so, it executes
    the PathogenFinder2 analysis using the pathogen_runner.py module. Otherwise, it passes
    the request to the next handler in the chain.
    """
    async def handle(self, analysis_name: str, db_folder: str, params: dict) -> asyncio.Task | None:
        """
        Handles PathogenFinder2 analysis requests.
        
        Args:
            analysis_name (str): The name of the analysis to perform.
            db_folder (str): The name of the database folder.
            params (dict): A dictionary of parameters specific to this analysis.
            
        Returns:
            asyncio.Task | None: A task for the running analysis, or None if not handled.
        """
        # Step 1: Check if this handler is responsible for the analysis.
        if analysis_name == "Pathogen_Finder2":
            # Step 2: If responsible, create and return a task for the specific workflow.
            return asyncio.create_task(self._run_pathogenfinder2_workflow(params))
        else:
            # Step 3: If not responsible, pass the request to the next handler in the chain.
            return await super().handle(analysis_name, db_folder, params)

    async def _run_pathogenfinder2_workflow(self, params: dict):
        """
        Runs the complete PathogenFinder2 analysis workflow.
        
        This function encapsulates the entire process for PathogenFinder2 analysis,
        including dependency checking, setup, execution, validation, and cleanup.
        
        Args:
            params (dict): A dictionary containing PathogenFinder2-specific parameters.
        """
        # Step 1: Announce the start of the workflow.
        self._context.logger.log_step("Pathogen_Finder2", "1_Start_Workflow", "PathogenFinder2 workflow initiated.")
        
        try:
            # Step 2: Setup - Check dependencies and prepare environment
            await self.setup()
            
            # Step 3: Execute the PathogenFinder2 analysis
            await self.execute()
            
            # Step 4: Validate the results
            await self.validate_results()
            
            # Step 5: Log completion
            self._context.logger.log_step("Pathogen_Finder2", "5_Workflow_Complete", "PathogenFinder2 workflow completed successfully.")
            
        except Exception as e:
            # Step 6: Handle any errors during the workflow
            self._context.logger.log_step("Pathogen_Finder2", "5_Workflow_Failed", f"PathogenFinder2 workflow failed: {str(e)}")
            raise

    async def setup(self):
        """
        Set up PathogenFinder2 environment and configuration.
        
        This method checks for required dependencies and sets up the configuration
        file with the correct input genome path.
        """
        # Step 1: Check for PathogenFinder2 dependencies
        self._context.logger.log_step("Pathogen_Finder2", "2_Check_Dependencies", "Checking PathogenFinder2 dependencies.")
        
        dependencies = ["prodigal", "python", "diamond"]
        missing_deps = []
        
        for dep in dependencies:
            if subprocess.run(["which", dep], capture_output=True, text=True).returncode != 0:
                missing_deps.append(dep)
        
        if missing_deps:
            error_msg = f"Missing PathogenFinder2 dependencies: { ', '.join(missing_deps)}"
            self._context.logger.log_step("Pathogen_Finder2", "2_Dependencies_Missing", error_msg)
            raise RuntimeError(error_msg)
        
        self._context.logger.log_step("Pathogen_Finder2", "2_Dependencies_OK", "All PathogenFinder2 dependencies found.")
        
        # Step 2: Set up configuration file
        self._context.logger.log_step("Pathogen_Finder2", "3_Setup_Config", "Setting up PathogenFinder2 configuration.")
        
        # Create output directory
        output_dir = self._context.results_dir / "Pathogen_Finder2"
        output_dir.mkdir(exist_ok=True)
        
        # Set up configuration file path
        config_file = output_dir / "config.json"
        
        # Create configuration with input genome path and required sections
        config_data = {
            "Misc Parameters": {
                "Notes": "This is a base config file",
                "Results Folder": "",
                "Name": "PathogenFinder2 Run",
                "Actions": ["inference"],
                "Report Results": "file",
                "Project Name": "PathogenFinder2",
                "Prodigal Path": "prodigal",
                "protT5 Path": "protT5",
                "protT5 Model": "Rostlab/ProstT5",
                "Diamond Path": "diamond"
            },
            "input_genome": str(self._context.genome_db_path),
            "output_dir": str(output_dir),
            "database_dir": str(Path.cwd() / "database" / "Pathogenfinder"),
            "Train Parameters": {
                "Optimizer": "NAdam",
                "Imbalance Sample": False,
                "Imbalance Weight": False,
                "Learning Rate": 1e-4,
                "Norm Scale": 1e-6,
                "Stochastic Depth Prob": 0.2,
                "Epochs": 5,
                "Lr Scheduler": "ReduceLROnPlateau",
                "Warm Up": 5,
                "Weight Decay": 1e-4,
                "Lr End": 1,
                "Mix Prec": True,
                "Asynchronity": True,
                "Num Workers": 8,
                "Bucketing": 12,
                "Stratified": True,
                "Data Sample": False,
                "Early Stopping": False,
                "Save Model": "best_epoch",
                "Prot Dim Split": False,
                "Loss Function": "BCEWithLogitsLoss",
                "Train DF": "/path/to/metadata_train.tsv",
                "Train Loc": "/path/to/folder_with_data/",
                "Validation DF": "/path/to/metadata_val.tsv",
                "Validation Loc": "/path/to/folder_with_data/",
                "Train Results": "dictionary",
                "Memory Report": False,
                "Wandb Report": False,
                "Results dir": "/path/to/folderoutput/"
            }
        }
        
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        self._context.logger.log_step("Pathogen_Finder2", "3_Config_Setup", f"Configuration file created at: {config_file}")
        
        # Add logging to validate configuration structure
        self._context.logger.log_step("Pathogen_Finder2", "3_Config_Debug",
                                     f"Config data structure: {json.dumps(config_data, indent=2)}")
        self._context.logger.log_step("Pathogen_Finder2", "3_Config_Debug",
                                     f"Config file exists: {config_file.exists()}")
        self._context.logger.log_step("Pathogen_Finder2", "3_Config_Debug",
                                     f"Config file content preview: {config_file.read_text()[:200]}...")
        
        # Store config file path for later use
        self.config_file = config_file
        self.output_dir = output_dir
        # Use the original genome file path, not the BLAST database path
        self.genome_file = Path.cwd() / "genome" / "GCF_000523395.1_10982.fasta_genomic.fna"

    async def execute(self):
        """
        Execute PathogenFinder2 analysis using pathogen_runner.py.
        
        This method runs the PathogenFinder2 analysis using the configuration
        file created during setup.
        """
        # Step 1: Announce execution start
        self._context.logger.log_step("Pathogen_Finder2", "4_Start_Execution", "Starting PathogenFinder2 execution.")
        
        # Step 2: Execute PathogenFinder2 using pathogen_runner
        # Note: This assumes pathogen_runner has a function to run PathogenFinder2
        # For now, we'll use a placeholder command structure
        # Use the inference configuration file we created
        inference_config = Path.cwd() / "database" / "Pathogenfinder" / "configs" / "config_inference.json"
        
        command = [
            "pathogenfinder2",  # This should be the actual PathogenFinder2 command
            "predict",  # Use the predict subcommand
            "-i", str(self.genome_file),  # Input genome file
            "-o", str(self.output_dir),  # Output directory
            "-f", "genome",  # Format is genome
            "-c", str(inference_config),  # Use the inference config file
            "--prodigalPath", "prodigal",  # Path to prodigal
            "--protT5Path", "protT5",  # Path to protT5
            "--diamondPath", "diamond"  # Path to diamond
        ]
        
        # Log the command being executed
        self._context.logger.log_step("Pathogen_Finder2", "4_Command_Debug",
                                     f"Executing command: {' '.join(command)}")
        self._context.logger.log_step("Pathogen_Finder2", "4_Command_Debug",
                                     f"Config file path: {self.config_file}")
        self._context.logger.log_step("Pathogen_Finder2", "4_Command_Debug",
                                     f"Config file exists: {self.config_file.exists()}")
        
        success, stdout, stderr = await pathogen_runner.run_command_async(command)
        
        # Log execution results
        self._context.logger.log_step("Pathogen_Finder2", "4_Command_Debug",
                                     f"Command execution success: {success}")
        self._context.logger.log_step("Pathogen_Finder2", "4_Command_Debug",
                                     f"Command stdout: {stdout[:500]}...")
        self._context.logger.log_step("Pathogen_Finder2", "4_Command_Debug",
                                     f"Command stderr: {stderr[:500]}...")
        
        if not success:
            error_msg = f"PathogenFinder2 execution failed: {stderr}"
            self._context.logger.log_step("Pathogen_Finder2", "4_Execution_Failed", error_msg)
            raise RuntimeError(error_msg)
        
        # Step 4: Log execution results
        self._context.logger.log_step("Pathogen_Finder2", "4_Execution_Complete",
                                     f"PathogenFinder2 execution completed. Output: {stdout}")
        
        # Store results for validation
        self.execution_results = stdout

    async def validate_results(self):
        """
        Validate PathogenFinder2 output.
        
        This method checks that the analysis produced valid output files
        and parses the results for storage.
        """
        # Step 1: Announce validation start
        self._context.logger.log_step("Pathogen_Finder2", "5_Start_Validation", "Starting PathogenFinder2 result validation.")
        
        # Step 2: Check for expected output files
        expected_files = [
            "pathogenfinder_results.tsv",
            "pathogenfinder_summary.txt"
        ]
        
        missing_files = []
        for filename in expected_files:
            filepath = self.output_dir / filename
            if not filepath.exists():
                missing_files.append(filename)
        
        if missing_files:
            error_msg = f"Missing expected output files: { ', '.join(missing_files)}"
            self._context.logger.log_step("Pathogen_Finder2", "5_Validation_Failed", error_msg)
            raise FileNotFoundError(error_msg)
        
        # Step 3: Parse and store results
        try:
            # Parse the main results file
            results_file = self.output_dir / "pathogenfinder_results.tsv"
            results_df = pd.read_csv(results_file, sep='\t')
            
            # Parse the summary file
            summary_file = self.output_dir / "pathogenfinder_summary.txt"
            with open(summary_file, 'r') as f:
                summary_content = f.read()
            
            # Store results in context
            self._context.results_data['pathogenfinder2'] = {
                'results': results_df.to_dict('records'),
                'summary': summary_content,
                'output_dir': str(self.output_dir)
            }
            
            # Step 4: Log validation success
            self._context.logger.log_step("Pathogen_Finder2", "5_Validation_Complete",
                                         f"PathogenFinder2 validation successful. Found {len(results_df)} results.")
            
        except Exception as e:
            error_msg = f"Failed to parse PathogenFinder2 results: {str(e)}"
            self._context.logger.log_step("Pathogen_Finder2", "5_Validation_Failed", error_msg)
            raise RuntimeError(error_msg)

    async def cleanup(self):
        """
        Clean up temporary files.
        
        This method removes temporary files created during the analysis.
        """
        # Step 1: Announce cleanup start
        self._context.logger.log_step("Pathogen_Finder2", "6_Start_Cleanup", "Starting PathogenFinder2 cleanup.")
        
        # Step 2: Remove temporary files
        temp_files_to_remove = [
            self.config_file,
        ]
        
        for temp_file in temp_files_to_remove:
            try:
                if temp_file.exists():
                    temp_file.unlink()
                    self._context.logger.log_step("Pathogen_Finder2", "6_File_Removed", f"Removed temporary file: {temp_file}")
            except Exception as e:
                self._context.logger.log_step("Pathogen_Finder2", "6_Cleanup_Warning",
                                             f"Warning: Could not remove {temp_file}: {str(e)}")
        
        # Step 3: Log cleanup completion
        self._context.logger.log_step("Pathogen_Finder2", "6_Cleanup_Complete", "PathogenFinder2 cleanup completed.")
