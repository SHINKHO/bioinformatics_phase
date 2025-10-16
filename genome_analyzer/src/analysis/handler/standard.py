
import asyncio
from pathlib import Path

from .base import AnalysisHandler
from analysis import blast_runner

class StandardAnalysisHandler(AnalysisHandler):
    """
    The default handler for all standard, single-step BLAST analyses.
    
    This handler is placed at the end of the chain. It processes any analysis
    that has not been handled by the preceding special-case handlers. It assumes
    the analysis is a straightforward BLAST search of a database against the genome.
    """
    async def handle(self, analysis_name: str, db_folder: str, params: dict) -> asyncio.Task | None:
        """
        Handles any analysis by treating it as a standard BLAST workflow.
        
        Since this is the last handler in the chain, it does not call `super().handle()`.
        It unconditionally creates a task to run the standard analysis workflow.
        
        Args:
            analysis_name (str): The name of the analysis (e.g., "Antimicrobial_Resistance").
            db_folder (str): The name of the database folder (e.g., "resfinder_db").
            params (dict): An empty dictionary (in this case), for interface compatibility.
            
        Returns:
            asyncio.Task: A task for the running standard analysis.
        """
        # Step 1: This handler is the default, so it always processes the request.
        # Create and return a task for the standard analysis workflow.
        return asyncio.create_task(self._run_other_analysis(db_folder, analysis_name))


    async def _run_other_analysis(self, db_folder: str, analysis_name: str):
        """
        Runs a standard BLAST-based analysis.
        
        This generic workflow combines all FASTA files from a given database folder,
        runs a single BLASTN search against the input genome, and saves the results.
        
        Related Functions:
        - blast_runner.run_blastn_async: Used to perform the BLAST search.
        
        Args:
            db_folder (str): The name of the database folder (e.g., "resfinder_db").
            analysis_name (str): The desired output name for the analysis.
        """
        # Step 1: Announce the start of the analysis.
        self._context.logger.log_step(analysis_name, "1_Start_Analysis", f"Starting {analysis_name} analysis.")
        
        # Step 2: Set up paths and find all database FASTA files.
        query_dir = Path.cwd() / "database" / db_folder
        output_dir = self._context.results_dir / analysis_name
        output_dir.mkdir(exist_ok=True)
        
        query_files = list(query_dir.rglob("*.f*a"))
        if not query_files:
            self._context.logger.log_step(analysis_name, "2_No_Fasta_Found", f"No FASTA files found in '{query_dir}', skipping.", extension="log")
            return

        # Step 3: Combine all found database files into a single query file.
        combined_query = output_dir / "combined_query.fasta"
        with open(combined_query, "w") as f_out:
            for f in query_files:
                f_out.write(f.read_text())
                
        # Step 4: Run the BLASTN search.
        output_path = output_dir / "blast_results.tsv"
        await blast_runner.run_blastn_async(combined_query, self._context.genome_db_path, output_path, {})
        
        # Step 5: Log the raw results and announce completion.
        with open(output_path, "r") as f:
            self._context.logger.log_step(analysis_name, "3_Blast_Results", f"BLAST search results for {analysis_name}:\n{f.read()}", extension="tsv")
        self._context.logger.log_step(analysis_name, "4_End_Analysis", f"Analysis '{analysis_name}' completed.")

