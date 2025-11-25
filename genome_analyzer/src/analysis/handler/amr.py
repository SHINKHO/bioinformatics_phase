import asyncio
import json
import pandas as pd
from pathlib import Path

from .base import AnalysisHandler
from analysis import blast_runner

class AMRHandler(AnalysisHandler):
    # ABRicate와 유사한 항균제 내성(AMR) 분석을 처리하는 핸들러.
    async def handle(self, analysis_name: str, db_folder: str, params: dict) -> asyncio.Task | None:
        # "Antimicrobial_Resistance" 분석 요청을 처리. In: analysis_name(str), db_folder(str), params(dict) / Out: asyncio.Task 또는 None
        if analysis_name == "Antimicrobial_Resistance":
            return asyncio.create_task(self._run_amr_workflow(db_folder, analysis_name))
        else:
            return await super().handle(analysis_name, db_folder, params)

    async def _run_amr_workflow(self, db_folder: str, analysis_name: str):
        # AMR 분석 워크플로우 실행. In: db_folder(str), analysis_name(str) / Out: None
        self._context.logger.log_step(analysis_name, "1_Start_AMR_Workflow", "AMR workflow initiated.")
        
        output_dir = self._context.results_dir / analysis_name
        output_dir.mkdir(exist_ok=True)

        query_dir = Path.cwd() / "database" / db_folder
        query_files = list(query_dir.rglob("*.f*a"))
        if not query_files:
            self._context.logger.log_step(analysis_name, "2_No_Fasta_Found", f"No FASTA files found in '{query_dir}', skipping.", extension="log")
            return

        combined_query = output_dir / "combined_query.fasta"
        with open(combined_query, "w") as f_out:
            for f in query_files:
                f_out.write(f.read_text())
                
        blast_results_path = output_dir / "blast_results.tsv"
        blast_options = {
            "perc_identity": 95, 
            "qcov_hsp_perc": 95, 
            "outfmt": "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp"
        }
        await blast_runner.run_blastn_async(combined_query, self._context.genome_db_path, blast_results_path, blast_options)
        
        with open(blast_results_path, "r") as f:
            self._context.logger.log_step(analysis_name, "3_Blast_Results", f"BLAST search results for {analysis_name}:\n{f.read()}", extension="tsv")

        try:
            df = pd.read_csv(blast_results_path, sep='\t', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovhsp'])
        except (pd.errors.EmptyDataError, FileNotFoundError):
            df = pd.DataFrame()

        summary_records = []
        if not df.empty:
            best_hits = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]

            for _, row in best_hits.iterrows():
                parts = row['qseqid'].split('_')
                gene = parts[0]
                database = db_folder
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

        self._context.results_data['amr'] = summary_records
        
        json_output_path = output_dir / "amr_summary.json"
        with open(json_output_path, "w") as f:
            json.dump(summary_records, f, indent=4)

        self._context.logger.log_step(analysis_name, "4_End_AMR_Workflow", f"AMR workflow completed. Found {len(summary_records)} genes.")
