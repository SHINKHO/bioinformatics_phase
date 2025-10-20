
from pathlib import Path

from .base import AnalysisHandler
from analysis import blast_runner

import pandas as pd
from Bio import SeqIO
import re

class MLSTHandler(AnalysisHandler):
    """
    A concrete handler for the special multi-step MLST workflow.
    """

    async def execute(self):
        self.logger.log_step(self.analysis_name, "1_Start_MLST_Workflow", "MLST workflow initiated.")

        loci_order = self._get_loci_order()
        if not loci_order:
            return # Error already logged in _get_loci_order

        housekeeping_blast_options = {"perc_identity": 90}
        allele_blast_options = {}

        extracted_genes_path = await self._extract_housekeeping_genes(loci_order, housekeeping_blast_options)
        if not extracted_genes_path:
            return

        profile, novel_alleles = await self._determine_allele_profile(extracted_genes_path, loci_order, allele_blast_options)

        st = self._find_sequence_type(profile, loci_order)

        results = {
            "filename": self.context.get('genome_id'),
            "scheme": self.context.get('species'),
            "st": st,
            "alleles": profile
        }
        self.context['results_data']['mlst'] = results
        
        self._write_output("results_json", results, write_type="json")
        if novel_alleles:
            novel_alleles_content = "".join([f">{rec.id} {rec.description}\n{rec.seq}\n" for rec in novel_alleles])
            self._write_output("novel_alleles", novel_alleles_content, write_type="text")

        self.logger.log_step(self.analysis_name, "6_End_MLST_Workflow", f"MLST workflow completed. Found ST: {st}, Profile: {profile}")

    def _get_loci_order(self) -> list:
        """Extracts loci order from the species profile file."""
        gene_dir = Path(self.inputs.get('species_db'))
        species = self.context.get('species')
        profile_file = gene_dir / f"{species}.txt"
        try:
            with open(profile_file, 'r') as f:
                header = f.readline().strip()
                return header.split('\t')[1:] # Skip ST column
        except (FileNotFoundError, IndexError) as e:
            self.logger.log_step(self.analysis_name, "Profile_Error", f"Could not read loci order from {profile_file}: {e}")
            return []

    async def _extract_housekeeping_genes(self, loci_order: list, blast_options: dict) -> Path | None:
        """Extracts housekeeping genes from the genome by blasting against known alleles."""
        gene_dir = Path(self.inputs.get('species_db'))
        temp_dir = Path(self.context.get('temp_dir'))
        probes_fasta = temp_dir / "mlst_probes.fasta"

        with open(probes_fasta, "w") as f_out:
            for locus in loci_order:
                record = next(SeqIO.parse(gene_dir / f"{locus}.tfa", "fasta"))
                SeqIO.write(record, f_out, "fasta")

        blast_result_path = temp_dir / "probes_vs_genome.tsv"
        await blast_runner.run_blastn_async(probes_fasta, self.context.get("genome_db_path"), blast_result_path, blast_options)
        self.logger.log_step(self.analysis_name, "3_Housekeeping_Gene_Blast_Results", (temp_dir / "probes_vs_genome.tsv").read_text(), extension="tsv")

        try:
            df = pd.read_csv(blast_result_path, sep='\t', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
            best_hits = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]
        except (pd.errors.EmptyDataError, KeyError):
            best_hits = pd.DataFrame()

        extracted_genes_path = temp_dir / "extracted_mlst_genes.fasta"
        with open(extracted_genes_path, "w") as f:
            if not best_hits.empty:
                for _, row in best_hits.iterrows():
                    locus = row['qseqid'].split('_')[0]
                    start, end = sorted((row['sstart'], row['send']))
                    strand = "plus" if row['sstart'] < row['send'] else "minus"
                    
                    success, stdout, stderr = await blast_runner.run_command_async(
                        ["blastdbcmd", "-db", str(self.context.get("genome_db_path")), "-entry", row['sseqid'], "-range", f"{start}-{end}", "-strand", strand]
                    )
                    if success and stdout and stdout.startswith('>'):
                        sequence = "".join(stdout.splitlines()[1:])
                        if sequence:
                            f.write(f">{locus}\n{sequence}\n")
                    else:
                        self.logger.log_step(self.analysis_name, f"Extraction_Failed_{locus}", f"blastdbcmd failed for {locus}.\nStderr: {stderr}")
        
        self.logger.log_step(self.analysis_name, "4_Extracted_Genes_Content", extracted_genes_path.read_text(), extension="fasta")
        return extracted_genes_path

    async def _determine_allele_profile(self, extracted_genes_path: Path, loci_order: list, blast_options: dict) -> tuple[dict, list]:
        """Determines the allele profile by blasting extracted genes against the allele database."""
        gene_dir = Path(self.inputs.get('species_db'))
        temp_dir = Path(self.context.get('temp_dir'))
        combined_alleles = temp_dir / "all_alleles.fasta"

        with open(combined_alleles, "w") as f_out:
            for locus_file in gene_dir.glob("*.tfa"):
                f_out.write(locus_file.read_text())
        
        allele_db_path = await blast_runner.create_blast_db_async(combined_alleles, temp_dir)
        
        blast_alleles_path = temp_dir / "genes_vs_alleles.tsv"
        await blast_runner.run_blastn_async(extracted_genes_path, allele_db_path, blast_alleles_path, blast_options)
        self.logger.log_step(self.analysis_name, "5_Allele_Blast_Results", blast_alleles_path.read_text(), extension="tsv")

        try:
            df_alleles = pd.read_csv(blast_alleles_path, sep='\t', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
            best_alleles = df_alleles.loc[df_alleles.groupby('qseqid')['bitscore'].idxmax()]
        except (pd.errors.EmptyDataError, KeyError):
            best_alleles = pd.DataFrame()

        profile = {}
        novel_alleles = []
        extracted_sequences = {rec.id: rec for rec in SeqIO.parse(extracted_genes_path, "fasta")}

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
                        novel_record.id = f"{locus}_novel_{self.context.get('genome_id')}"
                        novel_record.description = f"Novel allele for {locus}"
                        novel_alleles.append(novel_record)
            else:
                profile[locus] = "?"
        
        if novel_alleles:
            self.logger.log_step(self.analysis_name, "Novel_Alleles_Found", f"Found {len(novel_alleles)} novel alleles.")

        return profile, novel_alleles

    def _find_sequence_type(self, profile: dict, loci_order: list) -> str:
        """Finds the sequence type (ST) by matching the allele profile against the profile file."""
        gene_dir = Path(self.inputs.get('species_db'))
        species = self.context.get('species')
        profile_file = gene_dir / f"{species}.txt"
        
        profile_df = pd.read_csv(profile_file, sep='\t').astype(str)
        st = "Novel Profile"
        
        # Ensure the profile has values for all loci in the correct order
        profile_values = [profile.get(locus, "?") for locus in loci_order]

        for _, row in profile_df.iterrows():
            row_values = list(row[loci_order])
            if row_values == profile_values:
                st = f"ST{row['ST']}"
                break
        return st
