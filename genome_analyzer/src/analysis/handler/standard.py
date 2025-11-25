
import asyncio
from pathlib import Path

from .base import AnalysisHandler
from analysis import blast_runner

class StandardAnalysisHandler(AnalysisHandler):
    # 표준 단일 단계 BLAST 분석을 처리하는 기본 핸들러.
    async def handle(self, analysis_name: str, db_folder: str, params: dict) -> asyncio.Task | None:
        # 표준 BLAST 워크플로우로 분석을 처리. In: analysis_name(str), db_folder(str), params(dict) / Out: asyncio.Task
        return asyncio.create_task(self._run_other_analysis(db_folder, analysis_name))


    async def _run_other_analysis(self, db_folder: str, analysis_name: str):
        # 표준 BLAST 기반 분석 실행. In: db_folder(str), analysis_name(str) / Out: None
        self._context.logger.log_step(analysis_name, "1_Start_Analysis", f"Starting {analysis_name} analysis.")
        
        query_dir = Path.cwd() / "database" / db_folder
        output_dir = self._context.results_dir / analysis_name
        output_dir.mkdir(exist_ok=True)
        
        query_files = list(query_dir.rglob("*.f*a"))
        if not query_files:
            self._context.logger.log_step(analysis_name, "2_No_Fasta_Found", f"No FASTA files found in '{query_dir}', skipping.", extension="log")
            return

        combined_query = output_dir / "combined_query.fasta"
        with open(combined_query, "w") as f_out:
            for f in query_files:
                f_out.write(f.read_text())
                
        output_path = output_dir / "blast_results.tsv"
        await blast_runner.run_blastn_async(combined_query, self._context.genome_db_path, output_path, {})
        
        with open(output_path, "r") as f:
            self._context.logger.log_step(analysis_name, "3_Blast_Results", f"BLAST search results for {analysis_name}:\n{f.read()}", extension="tsv")
        self._context.logger.log_step(analysis_name, "4_End_Analysis", f"Analysis '{analysis_name}' completed.")
