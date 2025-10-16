
import asyncio
import re
from typing import List
import pandas as pd
from Bio import SeqIO

from .base import AnalysisHandler
from analysis import blast_runner
from config import DATABASE_ROOT

import json

class MLSTHandler(AnalysisHandler):
    """
    A concrete handler for the special multi-step MLST workflow.
    
    This handler checks if the requested analysis is "MLST". If so, it executes
    the complex MLST logic. Otherwise, it passes the request to the next handler.
    """
    async def handle(self, analysis_name: str, db_folder: str, params: dict) -> asyncio.Task | None:
        # Step 1: Check if this handler is responsible for the analysis.
        if analysis_name == "MLST":
            # Step 2: If responsible, create and return a task for the specific workflow.
            # Use species from context instead of parameters
            mlst_params = {
                'species': self._context.species,
                'gene_dir': DATABASE_ROOT / "MLST_DB" / self._context.species,
                'profile_file': DATABASE_ROOT / "MLST_DB" / self._context.species / f"{self._context.species}.txt",
                'loci_order': self._get_loci_order_from_profile(),
                'genome_id': self._context.genome_id
            }
            return asyncio.create_task(self._run_mlst_workflow(mlst_params))
        else:
            # Step 3: If not responsible, pass the request to the next handler in the chain.
            return await super().handle(analysis_name, db_folder, params)

    def _get_loci_order_from_profile(self) -> List[str]:
        """
        Extract loci order from the species profile file.
        
        Returns:
            List[str]: List of loci in the order they appear in the profile file.
        """
        profile_file = DATABASE_ROOT / "MLST_DB" / self._context.species / f"{self._context.species}.txt"
        try:
            with open(profile_file, 'r') as f:
                # Header line contains the loci order (first line)
                lines = f.readlines()
                if len(lines) > 0:
                    loci_line = lines[0].strip()
                    return loci_line.split('\t')[1:]  # Skip first column (ST)
                else:
                    self._context.logger.log_step("MLST", "Profile_Error", "Profile file is empty.")
                    return []
        except FileNotFoundError:
            self._context.logger.log_step("MLST", "Profile_Error", f"Profile file not found: {profile_file}")
            return []
        except Exception as e:
            self._context.logger.log_step("MLST", "Profile_Error", f"Error reading profile file: {str(e)}")
            return []

    async def _run_mlst_workflow(self, mlst_params: dict):
        self._context.logger.log_step("MLST", "1_Start_MLST_Workflow", "MLST workflow initiated.")
        output_dir = self._context.results_dir / "MLST"
        output_dir.mkdir(exist_ok=True)

        gene_dir = mlst_params['gene_dir']
        loci_order = mlst_params['loci_order']
        profile_file = mlst_params['profile_file']
        
        housekeeping_blast_options = {"perc_identity": 90}
        allele_blast_options = {} # Use defaults from blast_runner

        extracted_genes_path = await self._extract_housekeeping_genes(gene_dir, loci_order, housekeeping_blast_options)
        
        profile = await self._determine_allele_profile(extracted_genes_path, gene_dir, loci_order, allele_blast_options)
        
        st = self._find_sequence_type(profile, profile_file)

        results = {
            "filename": self._context.genome_id,
            "scheme": self._context.species,
            "st": st,
            "alleles": profile
        }
        self._context.results_data['mlst'] = results
        
        json_output_path = output_dir / "mlst_results.json"
        with open(json_output_path, "w") as f:
            json.dump(results, f, indent=4)

        self._context.logger.log_step("MLST", "6_End_MLST_Workflow", f"MLST workflow completed. Found ST: {st}, Profile: {profile}")

    async def _extract_housekeeping_genes(self, gene_dir, loci_order, blast_options):
        probes_fasta = self._context.temp_dir / "mlst_probes.fasta"
        with open(probes_fasta, "w") as f_out:
            for locus in loci_order:
                record = next(SeqIO.parse(gene_dir / f"{locus}.tfa", "fasta"))
                SeqIO.write(record, f_out, "fasta")

        blast_result_path = self._context.temp_dir / "probes_vs_genome.tsv"
        await blast_runner.run_blastn_async(probes_fasta, self._context.genome_db_path, blast_result_path, blast_options)
        
        with open(blast_result_path, "r") as f:
            self._context.logger.log_step("MLST","3_Housekeeping_Gene_Blast_Results", f"BLAST search results for housekeeping genes:\n{f.read()}", extension="tsv")

        try:
            df = pd.read_csv(blast_result_path, sep='\t', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
            best_hits = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]
        except (pd.errors.EmptyDataError, KeyError):
            best_hits = pd.DataFrame()

        extracted_genes_path = self._context.temp_dir / "extracted_mlst_genes.fasta"
        with open(extracted_genes_path, "w") as f:
            if not best_hits.empty:
                for _, row in best_hits.iterrows():
                    locus = row['qseqid'].split('_')[0]
                    start, end = sorted((row['sstart'], row['send']))
                    strand = "plus" if row['sstart'] < row['send'] else "minus"
                    
                    success, stdout, stderr = await blast_runner.run_command_async(
                        ["blastdbcmd", "-db", str(self._context.genome_db_path), "-entry", row['sseqid'], "-range", f"{start}-{end}", "-strand", strand]
                    )
                    if success and stdout and stdout.startswith('>'):
                        sequence = "".join(stdout.splitlines()[1:])
                        if sequence:
                            f.write(f">{locus}\n{sequence}\n")
                    else:
                        self._context.logger.log_step("MLST", f"Extraction_Failed_{locus}", f"blastdbcmd failed for {locus}.\nStderr: {stderr}")
        
        with open(extracted_genes_path, "r") as f:
            self._context.logger.log_step("MLST", "4_Extracted_Genes_Content", f"Content of extracted_mlst_genes.fasta:\n\n{f.read()}", extension="fasta")
        
        return extracted_genes_path

    async def _determine_allele_profile(self, extracted_genes_path, gene_dir, loci_order, blast_options):
        combined_alleles = self._context.temp_dir / "all_alleles.fasta"
        with open(combined_alleles, "w") as f_out:
            for locus_file in gene_dir.glob("*.tfa"):
                f_out.write(locus_file.read_text())
        
        allele_db_path = await blast_runner.create_blast_db_async(combined_alleles, self._context.temp_dir)
        
        blast_alleles_path = self._context.temp_dir / "genes_vs_alleles.tsv"
        await blast_runner.run_blastn_async(extracted_genes_path, allele_db_path, blast_alleles_path, blast_options)
        
        with open(blast_alleles_path, "r") as f:
            self._context.logger.log_step("MLST", "5_Allele_Blast_Results", f"BLAST results for allele determination:\n{f.read()}", extension="tsv")

        try:
            df_alleles = pd.read_csv(blast_alleles_path, sep='\t', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
            best_alleles = df_alleles.loc[df_alleles.groupby('qseqid')['bitscore'].idxmax()]
        except (pd.errors.EmptyDataError, KeyError):
            best_alleles = pd.DataFrame()

        profile = {}
        novel_alleles = []
        output_dir = self._context.results_dir / "MLST"

        extracted_sequences = {record.id: record for record in SeqIO.parse(extracted_genes_path, "fasta")}

        for locus in loci_order:
            hit = best_alleles[best_alleles['qseqid'] == locus] if not best_alleles.empty else pd.DataFrame()
            if not hit.empty:
                pident = hit.iloc[0]['pident']
                if pident >= 100.0:
                    allele_id = hit.iloc[0]['sseqid']
                    allele_num_match = re.search(r'(\d+)', allele_id)
                    profile[locus] = allele_num_match.group(1) if allele_num_match else "?"
                else:
                    profile[locus] = f"novel({pident:.2f}%)"
                    if locus in extracted_sequences:
                        novel_record = extracted_sequences[locus]
                        novel_record.id = f"{locus}_novel"
                        novel_record.description = f"Novel allele for {locus}"
                        novel_alleles.append(novel_record)
            else:
                profile[locus] = "?"

        if novel_alleles:
            novel_alleles_path = output_dir / "novel_alleles.fasta"
            with open(novel_alleles_path, "w") as f:
                SeqIO.write(novel_alleles, f, "fasta")
            self._context.logger.log_step("MLST", "Novel_Alleles_Found", f"Found {len(novel_alleles)} novel alleles. Saved to {novel_alleles_path}")

        return profile

    def _find_sequence_type(self, profile, profile_file):
        profile_df = pd.read_csv(profile_file, sep='\t').astype(str)
        st = "Novel Profile"
        for _, row in profile_df.iterrows():
            if list(row[profile_df.columns[1:]]) == [profile.get(locus, "?") for locus in profile_df.columns[1:]]:
                st = f"ST{row['ST']}"
                break
        return st
