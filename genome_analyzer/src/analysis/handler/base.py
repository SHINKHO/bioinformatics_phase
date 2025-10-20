
import asyncio
from abc import ABC, abstractmethod
from typing import Dict
import json
from pathlib import Path

from analysis.context import WorkflowContext
from analysis import blast_runner

# --- Handler ABC ---

class AnalysisHandler(ABC):
    """
    Abstract Base Class for all analysis handlers.
    
    This class provides a standardized structure for all analysis steps. It
    handles context management, I/O path resolution, and provides abstract
    methods for file reading and writing, ensuring consistency across the pipeline.
    """
    def __init__(self, context: WorkflowContext, step_config: dict):
        """
        Initializes the handler, resolving paths and creating directories.
        
        Args:
            context (WorkflowContext): The shared workflow context.
            step_config (dict): The configuration dictionary for this specific step.
        """
        self.context = context
        self.step_config = step_config
        self.inputs = {k: self.context.resolve_path(v) for k, v in step_config.get("inputs", {}).items()}
        self.outputs = {k: self.context.resolve_path(v) for k, v in step_config.get("outputs", {}).items()}
        self.parameters = step_config.get("parameters", {})
        self.logger = self.context.get("logger")
        self.analysis_name = self.parameters.get("analysis_name", self.step_config.get("name", "UnnamedAnalysis"))

        # Ensure all parent directories for output files exist.
        for path in self.outputs.values():
            Path(path).parent.mkdir(parents=True, exist_ok=True)

    def _read_input(self, input_name: str, read_type: str = 'text'):
        """
        Reads data from a declared input file.

        This is a common tool to abstract file reading. It can be extended
        to handle various file types like 'json', 'fasta', etc.

        Args:
            input_name (str): The key of the input file as defined in pipeline.json.
            read_type (str, optional): The type of file to read ('text', 'json'). Defaults to 'text'.

        Returns:
            The content of the file, or None if not found.
        """
        file_path = self.inputs.get(input_name)
        if not file_path or not Path(file_path).exists():
            self.logger.log_step(self.analysis_name, f"InputError", f"Input '{input_name}' not found at {file_path}")
            return None
        
        if read_type == 'json':
            with open(file_path, 'r') as f:
                return json.load(f)
        # Add other read types as needed (e.g., 'fasta', 'binary')
        else: # Default to text
            with open(file_path, 'r') as f:
                return f.read()

    def _write_output(self, output_name: str, data, write_type: str = 'text'):
        """
        Writes data to a declared output file.

        This is a common tool to abstract file writing. It can be extended
        to handle various data types and formats.

        Args:
            output_name (str): The key of the output file as defined in pipeline.json.
            data: The data to write to the file.
            write_type (str, optional): The format to write the data in ('text', 'json'). Defaults to 'text'.
        """
        file_path = self.outputs.get(output_name)
        if not file_path:
            self.logger.log_step(self.analysis_name, f"OutputError", f"Output '{output_name}' not defined in pipeline config.")
            return

        if write_type == 'json':
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=4)
        # Add other write types as needed
        else: # Default to text
            with open(file_path, 'w') as f:
                f.write(data)
        self.logger.log_step(self.analysis_name, f"OutputWritten", f"Output '{output_name}' saved to {file_path}")

    @abstractmethod
    async def execute(self):
        """
        The core analysis logic to be implemented by subclasses.
        
        This method will be called by the PipelineService and should contain
        the primary logic for the analysis step.
        """
        pass

    async def _run_blast(self, query_file: Path, output_file: Path, blast_options: dict = None):
        """
        Runs a BLASTn search and logs the results.

        Args:
            query_file (Path): The path to the query FASTA file.
            output_file (Path): The path to the BLAST output file.
            blast_options (dict, optional): BLAST options. Defaults to None.
        """
        if blast_options is None:
            blast_options = {}
            
        await blast_runner.run_blastn_async(query_file, self.context.get("genome_db_path"), output_file, blast_options)
        
        with open(output_file, "r") as f:
            self.logger.log_step(self.analysis_name, "3_Blast_Results", f"BLAST search results for {self.analysis_name}:\n{f.read()}", extension="tsv")
