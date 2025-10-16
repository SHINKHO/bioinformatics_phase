
import asyncio
from abc import ABC, abstractmethod
from pathlib import Path
from dataclasses import dataclass

from logger import Logger

# --- Data Context ---

@dataclass
class AnalysisContext:
    """
    A data class to hold shared data, tools, and configurations needed by all handlers.
    This avoids passing many individual arguments through the handler chain.
    
    Attributes:
        genome_db_path (Path): The file path to the BLAST database of the input genome.
        results_dir (Path): The main directory where final results are stored.
        temp_dir (Path): A directory for intermediate files.
        logger (Logger): The instance of the logger for detailed step-logging.
        verbose (bool): Flag to enable verbose console output.
        results_data (dict): A dictionary to store the final results from all analyses.
        genome_id (str): The genome identifier extracted from folder structure.
        species (str): The species name extracted from folder structure.
    """
    genome_db_path: Path
    results_dir: Path
    temp_dir: Path
    logger: Logger
    verbose: bool
    results_data: dict
    genome_id: str
    species: str

# --- Handler ABC ---

class AnalysisHandler(ABC):
    """
    Abstract Base Class for all analysis handlers.
    
    This class defines the common interface for all handlers in the chain of
    responsibility. It includes methods to link handlers together (`set_next`)
    and to process an analysis request (`handle`).
    
    Attributes:
        _next_handler (AnalysisHandler | None): The next handler in the chain.
        _context (AnalysisContext): Shared data and tools.
    """
    def __init__(self, context: AnalysisContext):
        self._next_handler: AnalysisHandler | None = None
        self._context = context

    def set_next(self, handler: 'AnalysisHandler') -> 'AnalysisHandler':
        """
        Links this handler to the next one in the chain.
        
        This allows chaining multiple handlers together. E.g., chain.set_next(handler1).set_next(handler2).
        
        Args:
            handler (AnalysisHandler): The next handler object to link to.
            
        Returns:
            AnalysisHandler: The next handler, to allow for fluent chaining.
        """
        self._next_handler = handler
        return handler

    @abstractmethod
    async def handle(self, analysis_name: str, db_folder: str, params: dict) -> asyncio.Task | None:
        """
        Handles an analysis request.
        
        If the handler can process this `analysis_name`, it does so and returns
        an asyncio.Task. Otherwise, it passes the request to the next handler
        in the chain.
        
        Args:
            analysis_name (str): The name of the analysis to perform (e.g., "MLST").
            db_folder (str): The name of the database folder for this analysis.
            params (dict): A dictionary of parameters specific to this analysis.
            
        Returns:
            asyncio.Task | None: A task for the running analysis, or None if not handled.
        """
        if self._next_handler:
            return await self._next_handler.handle(analysis_name, db_folder, params)
        return None
