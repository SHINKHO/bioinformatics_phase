
import asyncio
from pathlib import Path

from .base import AnalysisHandler
from analysis import utils

class StandardAnalysisHandler(AnalysisHandler):
    """
    The default handler for all standard, single-step BLAST analyses.
    
    This handler is placed at the end of the chain. It processes any analysis
    that has not been handled by the preceding special-case handlers. It assumes
    the analysis is a straightforward BLAST search of a database against the genome.
    """
    async def execute(self):
        self.logger.log_step(self.analysis_name, "1_Start_Analysis", f"Starting {self.analysis_name} analysis.")
        
        query_dir = Path(self.inputs.get("db_folder"))
        output_dir = Path(self.outputs.get("blast_results")).parent
        
        combined_query = output_dir / "combined_query.fasta"
        if not utils.combine_fasta_files(query_dir, combined_query, self.logger, self.analysis_name):
            return

        output_path = Path(self.outputs.get("blast_results"))
        blast_options = self.parameters.get("blast_options", {})
        await self._run_blast(combined_query, output_path, blast_options)
        
        self.logger.log_step(self.analysis_name, "4_End_Analysis", f"Analysis '{self.analysis_name}' completed.")

