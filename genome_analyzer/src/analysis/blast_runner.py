# src/analysis/blast_runner.py (최종 수정 코드)
import asyncio
from pathlib import Path
from typing import Dict, Tuple

async def run_command_async(command: list) -> Tuple[bool, str, str]:
    """
    Runs a command asynchronously using asyncio.create_subprocess_exec.
    This is generally more robust for handling arguments than create_subprocess_shell.
    """
    # Ensure all command parts are strings
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
        # This happens if the command itself (e.g., 'blastn') isn't found
        tool = command[0]
        return False, "", f"Error: Command '{tool}' not found. Is it installed and in your PATH?"
    except Exception as e:
        return False, "", f"An unexpected error occurred: {e}"


async def create_blast_db_async(fasta_file: Path, db_output_dir: Path) -> Path:
    """Creates a BLAST database from a FASTA file if it doesn't already exist."""
    db_name = db_output_dir / fasta_file.stem
    
    # Check if DB files already exist to avoid re-creation
    if not any(db_name.with_suffix(s).exists() for s in ['.nin', '.nhr', '.nsq']):
        command = [
            "makeblastdb",
            "-in", str(fasta_file),
            "-dbtype", "nucl",
            "-out", str(db_name),
            "-parse_seqids" # Important for FASTA files with complex headers
        ]
        success, stdout, stderr = await run_command_async(command)
        if not success:
            raise RuntimeError(f"makeblastdb failed: {stderr}")
            
    return db_name


async def run_blastn_async(query_file: Path, db_path: Path, output_file: Path, options: Dict):
    """Runs a BLASTN search with specified options."""
    # Default options
    default_opts = {
        "outfmt": "6",
        "perc_identity": 95,
        "qcov_hsp_perc": 95
    }
    # User-provided options override defaults
    final_opts = {**default_opts, **options}

    command = [
        "blastn",
        "-query", str(query_file),
        "-db", str(db_path),
        "-out", str(output_file)
    ]
    # Add options to the command list
    for key, value in final_opts.items():
        command.extend([f"-{key}", str(value)])

    success, stdout, stderr = await run_command_async(command)
    if not success:
        # If BLAST fails but produces some output, we might not want to raise an error
        # But if stderr has content, it's likely a real problem.
        if stderr:
            raise RuntimeError(f"blastn failed: {stderr}")