import asyncio
import json
import pandas as pd
from pathlib import Path

from .base import AnalysisHandler
from analysis import utils

class AMRHandler(AnalysisHandler):
    """
    A concrete handler for Antimicrobial Resistance (AMR) analysis, mimicking ABRicate.
    """
    async def execute(self):
        self.logger.log_step(self.analysis_name, "1_Start_AMR_Workflow", "AMR workflow initiated.")
        
        output_dir = Path(self.outputs.get("summary_json")).parent
        output_dir.mkdir(exist_ok=True) # Ensure output dir exists

        query_dir = Path(self.inputs.get("amr_db_folder"))
        
        combined_query = output_dir / "combined_query.fasta"
        if not utils.combine_fasta_files(query_dir, combined_query, self.logger, self.analysis_name):
            return
                
        blast_results_path = output_dir / "blast_results.tsv" # This is an intermediate file now
        blast_options = {
            "perc_identity": 95, 
            "qcov_hsp_perc": 95, 
            "outfmt": "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp"
        }
        await self._run_blast(combined_query, blast_results_path, blast_options)

        try:
            df = pd.read_csv(blast_results_path, sep='\t', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovhsp'])
        except (pd.errors.EmptyDataError, FileNotFoundError):
            df = pd.DataFrame()

        summary_records = []
        if not df.empty:
            # ... (parsing logic remains the same)
            best_hits = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]

            for _, row in best_hits.iterrows():
                parts = row['qseqid'].split('_')
                gene = parts[0]
                database = self.inputs.get("amr_db_folder")
                accession = parts[2] if len(parts) > 2 else 'N/A'
                product = "N/A"

                summary_records.append({
                    "GENE": gene,
                    "SEQUENCE": row['sseqid'],
                    "%COVERAGE": row['qcovhsp'],
                    "%IDENTITY": row['pident'],
                    "DATABASE": database,
                    "ACCESSION": accession,
                    "PRODUCT": product
                })

        # The main output is the summary JSON
        self.context['results_data']['amr'] = summary_records # Store in context
        
        self._write_output("summary_json", summary_records, write_type="json")

        self.logger.log_step(self.analysis_name, "4_End_AMR_Workflow", f"AMR workflow completed. Found {len(summary_records)} genes.")
