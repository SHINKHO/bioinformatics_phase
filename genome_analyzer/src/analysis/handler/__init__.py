"""
This package contains the various analysis handlers for the genomics pipeline.

Each handler is responsible for a specific type of analysis, following the
Chain of Responsibility pattern. The handlers are dynamically loaded and
chained together by the AnalysisManager.

- `base.py`: Defines the abstract `AnalysisHandler` and the `AnalysisContext` data class.
- `mlst.py`: Implements the multi-step MLST (Multi-Locus Sequence Typing) workflow.
- `pathogen_finder.py`: Implements the PathogenFinder2 workflow.
- `standard.py`: A generic handler for standard single-step BLAST analyses.
"""
from .base import AnalysisContext, AnalysisHandler
from .mlst import MLSTHandler
from .pathogen_finder import PathogenFinder2Handler
from .standard import StandardAnalysisHandler
from .amr import AMRHandler

from .amr import AMRHandler

__all__ = [
    "AnalysisContext",
    "AnalysisHandler",
    "MLSTHandler",
    "PathogenFinder2Handler",
    "StandardAnalysisHandler",
    "AMRHandler",
]