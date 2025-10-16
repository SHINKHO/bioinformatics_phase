"""
Asynchronous Pathogen Finder2 Command Wrappers

This module provides high-level asynchronous wrappers for common genomicepidemiology/PathogenFinder2
command-line tools. It uses `asyncio.create_subprocess_exec` for robust and
efficient execution of external bioinformatics tools like `blastn` and `makeblastdb`.

This approach allows multiple PathogenFinder setup to process running concurrently, significantly
speeding up the pipeline.
"""
import asyncio
from pathlib import Path
from typing import Dict, Tuple

async def run_command_async(command: list) -> Tuple[bool, str, str]:
    """
    Runs a given command in a subprocess asynchronously.

    This function serves as a general-purpose asynchronous command runner. It captures
    stdout and stderr, and returns them along with a success status. It's designed
    to be more robust than `asyncio.create_subprocess_shell` by taking arguments
    as a list, which avoids shell injection issues.

    Args:
        command (list): A list of strings representing the command and its arguments
                        (e.g., ["blastn", "-query", "q.fasta", ...]).

    Returns:
        Tuple[bool, str, str]: A tuple containing:
            - bool: True if the command executed successfully (return code 0), False otherwise.
            - str: The content of stdout from the command.
            - str: The content of stderr from the command.
    """
    # Step 1: Ensure all parts of the command are strings for compatibility.
    cmd_str_list = [str(item) for item in command]

    try:
        # Step 2: Create the asynchronous subprocess.
        # stdout and stderr are piped to be captured.
        proc = await asyncio.create_subprocess_exec(
            *cmd_str_list,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        
        # Step 3: Wait for the process to complete and communicate results.
        stdout_bytes, stderr_bytes = await proc.communicate()
        
        # Step 4: Decode the byte output to strings.
        stdout = stdout_bytes.decode('utf-8', errors='ignore')
        stderr = stderr_bytes.decode('utf-8', errors='ignore')

        # Step 5: Return the success status and the captured output.
        return proc.returncode == 0, stdout, stderr

    except FileNotFoundError:
        # Step 6a: Handle error if the command executable itself is not found.
        tool = command[0]
        return False, "", f"Error: Command '{tool}' not found. Is it installed and in your PATH?"
    except Exception as e:
        # Step 6b: Handle any other unexpected errors during subprocess execution.
        return False, "", f"An unexpected error occurred: {e}"


async def create_blast_db_async(fasta_file: Path, db_output_dir: Path) -> Path:
    """
    Creates a BLAST database from a given FASTA file if it doesn't already exist.

    This function checks for the existence of BLAST database files (.nin, .nhr, .nsq)
    to avoid redundant database creation. If they don't exist, it calls `makeblastdb`.

    Related Functions:
    - run_command_async: Used to execute the `makeblastdb` command.

    Args:
        fasta_file (Path): The path to the input FASTA file.
        db_output_dir (Path): The directory where the BLAST database files will be stored.

    Returns:
        Path: The base path to the newly created BLAST database (without extension).

    Raises:
        RuntimeError: If the `makeblastdb` command fails.
    """
    # Step 1: Define the base name for the database from the input FASTA file.
    db_name = db_output_dir / fasta_file.stem
    
    # Step 2: Check if database files already exist to prevent re-creation.
    if not any(db_name.with_suffix(s).exists() for s in ['.nin', '.nhr', '.nsq']):
        # Step 3: If DB doesn't exist, construct the `makeblastdb` command.
        command = [
            "makeblastdb",
            "-in", str(fasta_file),
            "-dbtype", "nucl",
            "-out", str(db_name),
            "-parse_seqids" # Important for FASTA files with complex headers
        ]
        # Step 4: Execute the command asynchronously.
        success, stdout, stderr = await run_command_async(command)
        if not success:
            # Step 5: If command fails, raise an error with the details from stderr.
            raise RuntimeError(f"makeblastdb failed: {stderr}")
            
    # Step 6: Return the path to the database.
    return db_name


async def run_blastn_async(query_file: Path, db_path: Path, output_file: Path, options: Dict):
    """
    Runs a BLASTN search with a given set of options.

    This function constructs and executes a `blastn` command. It uses a set of
    default options for identity and coverage, which can be overridden by the
    `options` dictionary.

    Related Functions:
    - run_command_async: Used to execute the `blastn` command.

    Args:
        query_file (Path): Path to the query FASTA file.
        db_path (Path): Path to the BLAST database to search against.
        output_file (Path): Path where the BLAST results will be saved.
        options (Dict): A dictionary of additional BLAST options (e.g., {"evalue": 1e-5}).

    Raises:
        RuntimeError: If the `blastn` command fails and produces an error message.
    """
    # Step 1: Define default BLAST options.
    default_opts = {
        "outfmt": "6", # Tabular output format
        "perc_identity": 95,
        "qcov_hsp_perc": 95
    }
    # Step 2: Merge user-provided options with defaults. User options take precedence.
    final_opts = {**default_opts, **options}

    # Step 3: Construct the base `blastn` command.
    command = [
        "blastn",
        "-query", str(query_file),
        "-db", str(db_path),
        "-out", str(output_file)
    ]
    
    # Step 4: Append the final options to the command list.
    for key, value in final_opts.items():
        command.extend([f"-{key}", str(value)])

    # Step 5: Execute the command asynchronously.
    success, stdout, stderr = await run_command_async(command)
    if not success:
        # Step 6: If the command fails, check for content in stderr and raise an error.
        # A non-zero exit code might occur even with no hits, but stderr content
        # usually indicates a genuine problem.
        if stderr:
            raise RuntimeError(f"blastn failed: {stderr}")
