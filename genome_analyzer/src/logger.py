# src/logger.py
import logging
from pathlib import Path
import time
from datetime import datetime
import os

class Logger:
    """A simple logger to save step-by-step debug information."""

    def __init__(self, log_dir: Path):
        self.log_dir = log_dir
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.log_counts = {}

    def log_step(self, analysis_type: str, step_name: str, content: str, extension: str = "log"):
        """
        Logs the content to a file with a structured name.
        
        Args:
            analysis_type (str): The type of analysis (e.g., 'MLST', 'AMR').
            step_name (str): A descriptive name for the step (e.g., '1_species_identification').
            content (str): The text content to write to the log file.
            extension (str): The file extension for the log (e.g., 'log', 'tsv', 'fasta').
        """
        try:
            # Sanitize step_name to be a valid filename
            safe_step_name = "".join(c for c in step_name if c.isalnum() or c in ('_', '-')).rstrip()
            
            # Handle duplicate log entries for the same step
            log_key = (analysis_type, safe_step_name)
            if log_key not in self.log_counts:
                self.log_counts[log_key] = 0
            self.log_counts[log_key] += 1
            
            # Construct filename
            date_str = datetime.now().strftime("%Y-%m-%d")
            count = self.log_counts[log_key]
            log_file_name = f"{date_str}_{analysis_type}_{safe_step_name}_{count}.{extension}"
            log_file = self.log_dir / log_file_name

            with open(log_file, "w", encoding="utf-8") as f:
                f.write(content)
        except Exception as e:
            print(f"Failed to write log for step '{step_name}'. Error: {e}")

def setup_run_logger(run_log_dir: Path):
    """Sets up the main run logger for console output."""
    logger = logging.getLogger('main_logger')
    logger.setLevel(logging.INFO)
    
    # Prevent adding multiple handlers if function is called more than once
    if not logger.handlers:
        # Console handler
        c_handler = logging.StreamHandler()
        c_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        c_handler.setFormatter(c_format)
        logger.addHandler(c_handler)

    return logger