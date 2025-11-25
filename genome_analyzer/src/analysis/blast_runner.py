# 비동기 BLAST+ 커맨드 래퍼 모듈
import asyncio
from pathlib import Path
from typing import Dict, Tuple

async def run_command_async(command: list) -> Tuple[bool, str, str]:
    # 비동기 커맨드 실행. In: command(list) / Out: (성공여부, stdout, stderr)
    cmd_str_list = [str(item) for item in command]
    try:
        proc = await asyncio.create_subprocess_exec(
            *cmd_str_list,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        stdout_bytes, stderr_bytes = await proc.communicate()
        stdout = stdout_bytes.decode('utf-8', errors='ignore')
        stderr = stderr_bytes.decode('utf-8', errors='ignore')
        return proc.returncode == 0, stdout, stderr
    except FileNotFoundError:
        tool = command[0]
        return False, "", f"Error: Command '{tool}' not found. Is it installed and in your PATH?"
    except Exception as e:
        return False, "", f"An unexpected error occurred: {e}"


async def create_blast_db_async(fasta_file: Path, db_output_dir: Path) -> Path:
    # BLAST DB 생성. In: fasta_file(Path), db_output_dir(Path) / Out: DB 경로(Path)
    db_name = db_output_dir / fasta_file.stem
    if not any(db_name.with_suffix(s).exists() for s in ['.nin', '.nhr', '.nsq']):
        command = [
            "makeblastdb",
            "-in", str(fasta_file),
            "-dbtype", "nucl",
            "-out", str(db_name),
            "-parse_seqids"
        ]
        success, stdout, stderr = await run_command_async(command)
        if not success:
            raise RuntimeError(f"makeblastdb failed: {stderr}")
    return db_name


async def run_blastn_async(query_file: Path, db_path: Path, output_file: Path, options: Dict):
    # BLASTN 실행. In: query_file(Path), db_path(Path), output_file(Path), options(Dict) / Out: None
    default_opts = {
        "outfmt": "6",
        "perc_identity": 95,
        "qcov_hsp_perc": 95
    }
    final_opts = {**default_opts, **options}
    command = [
        "blastn",
        "-query", str(query_file),
        "-db", str(db_path),
        "-out", str(output_file)
    ]
    for key, value in final_opts.items():
        command.extend([f"-{key}", str(value)])
    success, stdout, stderr = await run_command_async(command)
    if not success:
        if stderr:
            raise RuntimeError(f"blastn failed: {stderr}")
