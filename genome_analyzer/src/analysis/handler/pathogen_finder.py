
import asyncio
import subprocess
from pathlib import Path
import pandas as pd
import json

from .base import AnalysisHandler
from analysis.context import WorkflowContext
from analysis import pathogen_runner

class PathogenFinder2Handler(AnalysisHandler):
    """
    A concrete handler for the PathogenFinder2 workflow.
    """
    async def execute(self):
        self.logger.log_step(self.analysis_name, "1_Start_Workflow", "PathogenFinder2 workflow initiated.")

        try:
            # --- Setup and Dependency Check ---
            self.logger.log_step(self.analysis_name, "2_Check_Dependencies", "Checking PathogenFinder2 dependencies.")
            dependencies = ["prodigal", "python", "diamond"]
            missing_deps = [dep for dep in dependencies if subprocess.run(["which", dep], capture_output=True).returncode != 0]
            if missing_deps:
                raise RuntimeError(f"Missing PathogenFinder2 dependencies: { ', '.join(missing_deps)}")
            self.logger.log_step(self.analysis_name, "2_Dependencies_OK", "All PathogenFinder2 dependencies found.")

            output_dir = Path(self.outputs.get("output_dir"))
            genome_file = Path(self.inputs.get("genome_file"))

            # --- Execute --- 
            self.logger.log_step(self.analysis_name, "4_Start_Execution", "Starting PathogenFinder2 execution.")
            inference_config = Path.cwd() / "database" / "Pathogenfinder" / "configs" / "config_inference.json"
            command = [
                "pathogenfinder2", "predict",
                "-i", str(genome_file),
                "-o", str(output_dir),
                "-f", "genome",
                "-c", str(inference_config),
                "--prodigalPath", "prodigal",
                "--protT5Path", "protT5", # Assuming these are in PATH
                "--diamondPath", "diamond"
            ]
            self.logger.log_step(self.analysis_name, "4_Command_Debug", f"Executing command: {' '.join(command)}")
            
            success, stdout, stderr = await pathogen_runner.run_command_async(command)

            if not success:
                raise RuntimeError(f"PathogenFinder2 execution failed: {stderr}")
            self.logger.log_step(self.analysis_name, "4_Execution_Complete", f"PathogenFinder2 execution completed.")

            # --- Validate Results ---
            self.logger.log_step(self.analysis_name, "5_Start_Validation", "Starting PathogenFinder2 result validation.")
            expected_files = ["pathogenfinder_results.tsv", "pathogenfinder_summary.txt"]
            missing_files = [f for f in expected_files if not (output_dir / f).exists()]
            if missing_files:
                raise FileNotFoundError(f"Missing expected output files: { ', '.join(missing_files)}")

            results_df = pd.read_csv(output_dir / "pathogenfinder_results.tsv", sep='\t')
            summary_content = (output_dir / "pathogenfinder_summary.txt").read_text()
            
            self.context['results_data']['pathogenfinder2'] = {
                'results': results_df.to_dict('records'),
                'summary': summary_content,
                'output_dir': str(output_dir)
            }
            self.logger.log_step(self.analysis_name, "5_Validation_Complete", f"PathogenFinder2 validation successful. Found {len(results_df)} results.")

        except Exception as e:
            self.logger.log_step(self.analysis_name, "5_Workflow_Failed", f"PathogenFinder2 workflow failed: {str(e)}")
            # Re-raise the exception to be caught by the pipeline service
            raise
