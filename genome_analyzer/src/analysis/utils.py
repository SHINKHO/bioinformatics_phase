# 분석 파이프라인을 지원하는 유틸리티 함수 모음.
import subprocess
from pathlib import Path
from Bio import SeqIO
from typing import Dict, Any, List

# Project-level module imports
from config import DATABASE_ROOT
from logger import Logger


def check_dependencies():
    # 필수 커맨드라인 도구들의 설치 여부를 확인. In: None / Out: None
    dependencies = ["blastn", "makeblastdb", "blastdbcmd", "prodigal", "python", "diamond"]
    for dep in dependencies:
        if subprocess.run(["which", dep], capture_output=True, text=True).returncode != 0:
            raise RuntimeError(
                f"Dependency '{dep}' not found in PATH. "
                "Please install NCBI BLAST+ and ensure it's in your system's PATH."
            )


def setup_mlst_parameters(genome_file: Path, logger: Logger) -> Dict[str, Any]:
    # 게놈 폴더 구조에서 종을 식별하고 MLST 분석을 위한 파라미터 준비. In: genome_file(Path), logger(Logger) / Out: MLST 파라미터(dict)
    try:
        parent_dir = genome_file.parent
        species_dir = parent_dir.name
        genome_id = parent_dir.parent.name
        
        logger.log_step("MLST", "1_1_Read_Folder_Structure",
                       f"Extracted from folder structure: genome_id='{genome_id}', species='{species_dir}'")
    except Exception as e:
        raise ValueError(
            f"Invalid genome folder structure. Expected format: {{genome_id}}/{{species}}/genome_file.fasta. "
            f"Got: '{genome_file}'. Error: {e}"
        )

    mlst_db_path = DATABASE_ROOT / "MLST_DB"
    potential_species = [d.name for d in mlst_db_path.iterdir() if d.is_dir()]
    
    if species_dir not in potential_species:
        raise ValueError(
            f"Species '{species_dir}' not found in MLST database. "
            f"Available species: {', '.join(potential_species)}"
        )

    species_db_dir = mlst_db_path / species_dir
    logger.log_step("MLST", "1_2_Species_Database_Found", f"MLST database found for species '{species_dir}' at '{species_db_dir}'")

    profile_file = next(species_db_dir.glob("*.txt"), None)
    if not profile_file:
        raise FileNotFoundError(f"MLST profile file (.txt) not found in '{species_db_dir}'")

    with open(profile_file, 'r') as f:
        header_line = f.readline().strip()
        loci_order = header_line.split('\t')[1:]

    return {
        "species": species_dir,
        "gene_dir": species_db_dir,
        "profile_file": profile_file,
        "loci_order": loci_order,
        "genome_id": genome_id
    }

def get_genome_name(genome_file: Path) -> str:
    # FASTA 파일에서 설명적인 이름을 추출. In: genome_file(Path) / Out: 첫 번째 레코드 ID(str)
    try:
        first_record = next(SeqIO.parse(genome_file, "fasta"))
        return first_record.id
    except StopIteration:
        raise ValueError(f"FASTA file '{genome_file}' is empty or not in a valid format.")