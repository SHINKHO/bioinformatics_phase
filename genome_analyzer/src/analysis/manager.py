# src/analysis/manager.py (최종 수정 코드)
import asyncio
import time
from pathlib import Path
import shutil
import pandas as pd
import re
from Bio import SeqIO
from datetime import datetime

from . import utils
from . import blast_runner
from reporting import reporter
from config import DATABASE_ROOT, ANALYSES_TO_RUN, BLAST_DB_DIR
from logger import Logger

class AnalysisManager:
    """Orchestrates the genome analysis workflow."""
    def __init__(self, genome_file: Path, results_dir: Path, verbose: bool = False):
        self.genome_file = genome_file
        self.results_dir = results_dir
        self.temp_dir = results_dir / "temp"
        self.verbose = verbose
        self.results_data = {}
        
        self.logs_dir = Path.cwd() / "logs"
        self.logger = Logger(self.logs_dir)

        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir.mkdir(exist_ok=True)

    def _log(self, message: str, level: str = "INFO"):
        """Prints a log message if in verbose mode."""
        if self.verbose:
            timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
            print(f"[{timestamp} - {level}] {message}")

    async def _run_mlst_workflow(self, genome_db_path: Path, mlst_params: dict):
        """Runs the complete MLST analysis sub-workflow."""
        self._log("Starting MLST workflow...")
        self.logger.log_step("MLST", "1_Start_MLST_Workflow", "MLST workflow initiated.")
        
        gene_dir = mlst_params['gene_dir']
        loci_order = mlst_params['loci_order']
        profile_file = mlst_params['profile_file']

        blast_options = {"perc_identity": 90}

        # 1. Extract MLST genes from genome
        self._log("Extracting MLST gene sequences from genome...")
        
        probes_fasta = self.temp_dir / "mlst_probes.fasta"
        with open(probes_fasta, "w") as f_out:
            for locus in loci_order:
                record = next(SeqIO.parse(gene_dir / f"{locus}.tfa", "fasta"))
                SeqIO.write(record, f_out, "fasta")

        blast_result_path = self.temp_dir / "probes_vs_genome.tsv"
        await blast_runner.run_blastn_async(probes_fasta, genome_db_path, blast_result_path, blast_options)
        
        with open(blast_result_path, "r") as f:
            blast_results = f.read()
            self.logger.log_step("MLST","3_Housekeeping_Gene_Blast_Results", f"BLAST search results for housekeeping genes:\n{blast_results}", extension="tsv")

        try:
            df = pd.read_csv(blast_result_path, sep='\t', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
            best_hits = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]
        except (pd.errors.EmptyDataError, KeyError):
            best_hits = pd.DataFrame()

        extracted_genes_path = self.temp_dir / "extracted_mlst_genes.fasta"
        with open(extracted_genes_path, "w") as f:
            if not best_hits.empty:
                for _, row in best_hits.iterrows():
                    locus = row['qseqid'].split('_')[0]
                    start, end = sorted((row['sstart'], row['send']))
                    strand = "plus" if row['sstart'] < row['send'] else "minus"
                    
                    # --- [최종 수정] blastdbcmd 실패에 대한 방어 코드 및 로깅 추가 ---
                    success, stdout, stderr = await blast_runner.run_command_async(
                        ["blastdbcmd", "-db", str(genome_db_path), "-entry", row['sseqid'], "-range", f"{start}-{end}", "-strand", strand]
                    )
                    # stdout이 비어있거나, FASTA 헤더(>)로 시작하지 않으면 실패로 간주
                    if success and stdout and stdout.startswith('>'):
                        sequence = "".join(stdout.splitlines()[1:])
                        if sequence: # 시퀀스 내용이 실제로 있는지 한번 더 확인
                            f.write(f">{locus}\n{sequence}\n")
                        else:
                             self._log(f"Extracted sequence for {locus} is empty. Skipping.", "WARN")
                             self.logger.log_step("MLST", f"Extraction_Warning_{locus}", f"Command output for {locus} was valid FASTA but sequence was empty.")
                    else:
                        self._log(f"Failed to extract sequence for {locus} using blastdbcmd. Skipping.", "WARN")
                        self.logger.log_step("MLST", f"Extraction_Failed_{locus}", f"blastdbcmd failed for {locus}.\nSuccess: {success}\nStderr: {stderr}\nStdout: {stdout}")
                    # --- 수정 끝 ---
        
        with open(extracted_genes_path, "r") as f:
            extracted_content = f.read()
            self.logger.log_step("MLST", "4_Extracted_Genes_Content", f"Content of extracted_mlst_genes.fasta:\n\n{extracted_content}", extension="fasta")

        # 2. Determine Allele types and ST
        self._log("Determining sequence type (ST)...")
        self.logger.log_step("MLST", "4_Start_ST_Determination", "Starting sequence type (ST) determination.")

        combined_alleles = self.temp_dir / "all_alleles.fasta"
        with open(combined_alleles, "w") as f_out:
            for locus_file in gene_dir.glob("*.tfa"):
                f_out.write(locus_file.read_text())
        
        allele_db_path = await blast_runner.create_blast_db_async(combined_alleles, self.temp_dir)
        blast_alleles_path = self.temp_dir / "genes_vs_alleles.tsv"
        
        await blast_runner.run_blastn_async(extracted_genes_path, allele_db_path, blast_alleles_path, blast_options)
        
        with open(blast_alleles_path, "r") as f:
            blast_results_alleles = f.read()
            self.logger.log_step("MLST", "5_Allele_Blast_Results", f"BLAST search results for allele determination:\n{blast_results_alleles}", extension="tsv")

        try:
            df_alleles = pd.read_csv(blast_alleles_path, sep='\t', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
        except pd.errors.EmptyDataError:
            df_alleles = pd.DataFrame() 

        if not df_alleles.empty:
            best_alleles = df_alleles.loc[df_alleles.groupby('qseqid')['bitscore'].idxmax()]
        else:
            best_alleles = pd.DataFrame()

        profile = {}
        for locus in loci_order:
            hit = best_alleles[best_alleles['qseqid'] == locus] if not best_alleles.empty else pd.DataFrame()
            if not hit.empty:
                allele_id = hit.iloc[0]['sseqid']
                allele_num_match = re.search(r'(\d+)', allele_id)
                if allele_num_match:
                    profile[locus] = allele_num_match.group(1)
                else:
                    profile[locus] = "?"
            else:
                profile[locus] = "?"

        profile_df = pd.read_csv(profile_file, sep='\t').astype(str)
        st = "Novel Profile"
        for _, row in profile_df.iterrows():
            if list(row[profile_df.columns[1:]]) == [profile.get(locus, "?") for locus in profile_df.columns[1:]]:
                st = f"ST{row['ST']}"
                break
        
        self.results_data['mlst'] = {"ST": st, "Profile": profile}
        self._log(f"MLST workflow completed. Found {st}.")
        self.logger.log_step("MLST", "6_End_MLST_Workflow", f"MLST workflow completed. Found ST: {st}, Profile: {profile}")

    async def _run_other_analysis(self, db_folder: str, analysis_name: str, genome_db_path: Path):
        self._log(f"Starting analysis: {analysis_name}...")
        self.logger.log_step(analysis_name, "1_Start_Analysis", f"Starting {analysis_name} analysis.")
        query_dir = DATABASE_ROOT / db_folder
        output_dir = self.results_dir / analysis_name
        output_dir.mkdir(exist_ok=True)
        query_files = list(query_dir.rglob("*.f*a"))
        if not query_files:
            self._log(f"No FASTA files found in '{query_dir}', skipping.", "WARN")
            self.logger.log_step(analysis_name, "2_No_Fasta_Found", f"No FASTA files found in '{query_dir}', skipping.", extension="log")
            return
        combined_query = output_dir / "combined_query.fasta"
        with open(combined_query, "w") as f_out:
            for f in query_files:
                f_out.write(f.read_text())
        output_path = output_dir / "blast_results.tsv"
        await blast_runner.run_blastn_async(combined_query, genome_db_path, output_path, {})
        with open(output_path, "r") as f:
            blast_results = f.read()
            self.logger.log_step(analysis_name, "3_Blast_Results", f"BLAST search results for {analysis_name}:\n{blast_results}", extension="tsv")
        self._log(f"Analysis '{analysis_name}' completed.")
        self.logger.log_step(analysis_name, "4_End_Analysis", f"Analysis '{analysis_name}' completed.")

    async def run_pipeline(self):
        start_time = time.time()
        print("===== Genome Analysis Pipeline Start =====")
        try:
            self._log("Step 1: Pre-flight checks and setup.")
            self.logger.log_step("Pipeline", "1_Pre-flight_Checks", "Starting pre-flight checks and setup.")
            utils.check_dependencies()
            self._log("Step 2: Creating BLAST database for the input genome.")
            self.logger.log_step("Pipeline", "3_Create_BLAST_DB", "Creating BLAST database for the input genome.")
            genome_db_path = await blast_runner.create_blast_db_async(self.genome_file, BLAST_DB_DIR)
            self._log(f"Genome BLAST DB created at '{genome_db_path}'.")
            self.logger.log_step("Pipeline", "4_BLAST_DB_Created", f"Genome BLAST DB created at '{genome_db_path}'.")
            self._log("Step 3: Running all analysis tasks concurrently.")
            self.logger.log_step("Pipeline", "5_Run_Concurrent_Analyses", "Running all analysis tasks concurrently.")
            mlst_params = utils.setup_mlst_parameters(self.genome_file, self.logger)
            self.results_data['mlst_params'] = mlst_params
            self._log(f"Species '{mlst_params['species']}' identified for MLST.")
            self.logger.log_step("Pipeline", "2_Species_Identification", f"Species '{mlst_params['species']}' identified.")
            tasks = []
            for db_folder, analysis_name in ANALYSES_TO_RUN.items():
                if analysis_name == "MLST":
                    tasks.append(self._run_mlst_workflow(genome_db_path, mlst_params))
                else:
                    tasks.append(self._run_other_analysis(db_folder, analysis_name, genome_db_path))
            await asyncio.gather(*tasks)
            self._log("All analysis tasks completed.")
            self.logger.log_step("Pipeline", "6_Concurrent_Analyses_Complete", "All analysis tasks completed.")
            self._log("Step 4: Generating final report.")
            self.logger.log_step("Pipeline", "7_Generate_Report", "Generating final report.")
            reporter.create_final_report(self.results_data, self.results_dir, self.genome_file)
            self._log("Final report generated.")
            self.logger.log_step("Pipeline", "8_Report_Generated", "Final report generated.")
        except (ValueError, FileNotFoundError, RuntimeError, Exception) as e:
            print(f"\n❌ PIPELINE FAILED: An error occurred.\n  -> {e}")
            self.logger.log_step("Pipeline", "9_Pipeline_Failed", f"PIPELINE FAILED: An error occurred.\n  -> {e}")
        finally:
            if self.temp_dir.exists():
                shutil.rmtree(self.temp_dir)
        end_time = time.time()
        print(f"\n🎉 Analysis complete in {end_time - start_time:.2f} seconds.")
        print(f"   Results are in '{self.results_dir}'")
        print("==========================================")