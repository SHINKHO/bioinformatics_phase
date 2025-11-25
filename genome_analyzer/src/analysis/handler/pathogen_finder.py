
import asyncio
import subprocess
from pathlib import Path
import pandas as pd
import json

from .base import AnalysisHandler
from analysis import pathogen_runner

class PathogenFinder2Handler(AnalysisHandler):
    # PathogenFinder2 워크플로우를 처리하는 핸들러.
    async def handle(self, analysis_name: str, db_folder: str, params: dict) -> asyncio.Task | None:
        # "Pathogen_Finder2" 분석 요청을 처리. In: analysis_name(str), db_folder(str), params(dict) / Out: asyncio.Task 또는 None
        if analysis_name == "Pathogen_Finder2":
            return asyncio.create_task(self._run_pathogenfinder2_workflow(params))
        else:
            return await super().handle(analysis_name, db_folder, params)

    async def _run_pathogenfinder2_workflow(self, params: dict):
        # PathogenFinder2 분석 워크플로우 실행. In: params(dict) / Out: None
        self._context.logger.log_step("Pathogen_Finder2", "1_Start_Workflow", "PathogenFinder2 workflow initiated.")
        
        try:
            await self.setup()
            await self.execute()
            await self.validate_results()
            self._context.logger.log_step("Pathogen_Finder2", "5_Workflow_Complete", "PathogenFinder2 workflow completed successfully.")
            
        except Exception as e:
            self._context.logger.log_step("Pathogen_Finder2", "5_Workflow_Failed", f"PathogenFinder2 workflow failed: {str(e)}")
            raise

    async def setup(self):
        # PathogenFinder2 환경 및 설정 구성. In: None / Out: None
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
        
        self._context.logger.log_step("Pathogen_Finder2", "3_Setup_Config", "Setting up PathogenFinder2 configuration.")
        
        output_dir = self._context.results_dir / "Pathogen_Finder2"
        output_dir.mkdir(exist_ok=True)
        
        config_file = output_dir / "config.json"
        
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
        
        self._context.logger.log_step("Pathogen_Finder2", "3_Config_Debug",
                                     f"Config data structure: {json.dumps(config_data, indent=2)}")
        self._context.logger.log_step("Pathogen_Finder2", "3_Config_Debug",
                                     f"Config file exists: {config_file.exists()}")
        self._context.logger.log_step("Pathogen_Finder2", "3_Config_Debug",
                                     f"Config file content preview: {config_file.read_text()[:200]}...")
        
        self.config_file = config_file
        self.output_dir = output_dir
        self.genome_file = Path.cwd() / "genome" / "GCF_000523395.1_10982.fasta_genomic.fna"

    async def execute(self):
        # pathogen_runner.py를 사용하여 PathogenFinder2 분석 실행. In: None / Out: None
        self._context.logger.log_step("Pathogen_Finder2", "4_Start_Execution", "Starting PathogenFinder2 execution.")
        
        inference_config = Path.cwd() / "database" / "Pathogenfinder" / "configs" / "config_inference.json"
        
        command = [
            "pathogenfinder2",
            "predict",
            "-i", str(self.genome_file),
            "-o", str(self.output_dir),
            "-f", "genome",
            "-c", str(inference_config),
            "--prodigalPath", "prodigal",
            "--protT5Path", "protT5",
            "--diamondPath", "diamond"
        ]
        
        self._context.logger.log_step("Pathogen_Finder2", "4_Command_Debug",
                                     f"Executing command: {' '.join(command)}")
        self._context.logger.log_step("Pathogen_Finder2", "4_Command_Debug",
                                     f"Config file path: {self.config_file}")
        self._context.logger.log_step("Pathogen_Finder2", "4_Command_Debug",
                                     f"Config file exists: {self.config_file.exists()}")
        
        success, stdout, stderr = await pathogen_runner.run_command_async(command)
        
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
        
        self._context.logger.log_step("Pathogen_Finder2", "4_Execution_Complete",
                                     f"PathogenFinder2 execution completed. Output: {stdout}")
        
        self.execution_results = stdout

    async def validate_results(self):
        # PathogenFinder2 결과 검증. In: None / Out: None
        self._context.logger.log_step("Pathogen_Finder2", "5_Start_Validation", "Starting PathogenFinder2 result validation.")
        
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
        
        try:
            results_file = self.output_dir / "pathogenfinder_results.tsv"
            results_df = pd.read_csv(results_file, sep='\t')
            
            summary_file = self.output_dir / "pathogenfinder_summary.txt"
            with open(summary_file, 'r') as f:
                summary_content = f.read()
            
            self._context.results_data['pathogenfinder2'] = {
                'results': results_df.to_dict('records'),
                'summary': summary_content,
                'output_dir': str(self.output_dir)
            }
            
            self._context.logger.log_step("Pathogen_Finder2", "5_Validation_Complete",
                                         f"PathogenFinder2 validation successful. Found {len(results_df)} results.")
            
        except Exception as e:
            error_msg = f"Failed to parse PathogenFinder2 results: {str(e)}"
            self._context.logger.log_step("Pathogen_Finder2", "5_Validation_Failed", error_msg)
            raise RuntimeError(error_msg)

    async def cleanup(self):
        # 임시 파일 정리. In: None / Out: None
        self._context.logger.log_step("Pathogen_Finder2", "6_Start_Cleanup", "Starting PathogenFinder2 cleanup.")
        
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
        
        self._context.logger.log_step("Pathogen_Finder2", "6_Cleanup_Complete", "PathogenFinder2 cleanup completed.")
