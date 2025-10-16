"""
Step-wise Debug Logger

This module provides a simple Logger class designed to save detailed, step-by-step
debugging information into organized log files. Each log entry is saved as a
separate file, named systematically based on the analysis type and step name.
"""
import logging
from pathlib import Path
from datetime import datetime

class Logger:
    """
    A simple logger to save step-by-step debug information into discrete files.
    
    This logger is used to trace the pipeline's execution, saving intermediate
    data, raw outputs, and command logs for debugging purposes.
    """

    def __init__(self, log_dir: Path):
        """
        Initializes the Logger.

        Args:
            log_dir (Path): The root directory where all log files will be stored.
        """
        # Step 1: Store the log directory path.
        self.log_dir = log_dir
        # Step 2: Create the log directory if it doesn't exist.
        self.log_dir.mkdir(parents=True, exist_ok=True)
        # Step 3: Initialize a dictionary to keep track of log counts for duplicate steps.
        self.log_counts = {}

    def log_step(self, analysis_type: str, step_name: str, content: str, extension: str = "log"):
        """
        Logs the given content to a file with a structured, unique name.
        
        Args:
            analysis_type (str): The type of analysis (e.g., 'MLST', 'AMR').
            step_name (str): A descriptive name for the step (e.g., '1_species_identification').
            content (str): The text content to write to the log file.
            extension (str): The file extension for the log (e.g., 'log', 'tsv', 'fasta').
        """
        try:
            # Step 1: Sanitize the step name to ensure it's a valid filename component.
            safe_step_name = "".join(c for c in step_name if c.isalnum() or c in ('_', '-')).rstrip()
            
            # Step 2: Handle potential duplicate log entries for the same step by appending a counter.
            log_key = (analysis_type, safe_step_name)
            if log_key not in self.log_counts:
                self.log_counts[log_key] = 0
            self.log_counts[log_key] += 1
            count = self.log_counts[log_key]
            
            # Step 3: Construct the full log filename with date, type, name, count, and extension.
            date_str = datetime.now().strftime("%Y-%m-%d")
            log_file_name = f"{date_str}_{analysis_type}_{safe_step_name}_{count}.{extension}"
            log_file = self.log_dir / log_file_name

            # Step 4: Write the content to the log file.
            with open(log_file, "w", encoding="utf-8") as f:
                f.write(content)
        except Exception as e:
            # Step 5: If logging fails, print an error to the console but do not halt the pipeline.
            print(f"Failed to write log for step '{step_name}'. Error: {e}")


def setup_run_logger(run_log_dir: Path):
    """
    (Deprecated/Unused) Sets up a main run logger for console output.
    
    Note: This function is currently not used in the main pipeline but is kept
    for potential future use. The current implementation uses a simpler, direct
    print statement controlled by the `--verbose` flag.

    Args:
        run_log_dir (Path): The directory where a main run log file could be stored.

    Returns:
        logging.Logger: A configured logger instance.
    """
    # Step 1: Get the logger instance.
    logger = logging.getLogger('main_logger')
    logger.setLevel(logging.INFO)
    
    # Step 2: Prevent adding duplicate handlers if the function is called multiple times.
    if not logger.handlers:
        # Step 3: Create and configure a console handler.
        c_handler = logging.StreamHandler()
        c_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        c_handler.setFormatter(c_format)
        logger.addHandler(c_handler)

    # Step 4: Return the configured logger.
    return logger
