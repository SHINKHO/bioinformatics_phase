"""
This package contains the various analysis handlers for the genomics pipeline.

Each handler is responsible for a specific type of analysis. The handlers are
dynamically loaded by the PipelineService based on the pipeline configuration.

- `base.py`: Defines the abstract `AnalysisHandler`.
- `mlst.py`: Implements the multi-step MLST (Multi-Locus Sequence Typing) workflow.
- `pathogen_finder.py`: Implements the PathogenFinder2 workflow.
- `standard.py`: A generic handler for standard single-step BLAST analyses.
"""
from .base import AnalysisHandler
from .mlst import MLSTHandler
from .pathogen_finder import PathogenFinder2Handler
from .standard import StandardAnalysisHandler
from .amr import AMRHandler

__all__ = [
    "AnalysisHandler",
    "MLSTHandler",
    "PathogenFinder2Handler",
    "StandardAnalysisHandler",
    "AMRHandler",
]