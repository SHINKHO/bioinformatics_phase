"""
Workflow Context for the Pipeline

This module defines the WorkflowContext class, which is a data container
that holds and manages state throughout the pipeline's execution.
"""
from __future__ import annotations
import re
from pathlib import Path
from typing import Any, Dict

class WorkflowContext:
    """
    A data class to hold shared data, tools, and configurations.
    
    This object is passed between pipeline steps to provide access to global
    configuration and the results of previous steps. It also provides
    helper methods for resolving path templates.
    """
    def __init__(self, initial_context: Dict[str, Any]):
        self._data = initial_context

    def __setitem__(self, key: str, value: Any):
        self._data[key] = value

    def __getitem__(self, key: str) -> Any:
        return self._data[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self._data.get(key, default)

    def resolve_path(self, path_template: str) -> str:
        """
        Resolves a path template string with values from the context.
        
        Example:
            resolve_path("{context.results_dir}/{context.species}/report.txt")
            -> "/path/to/results/klebsiella/report.txt"
        """
        # Regex to find all placeholders like {context.key}
        placeholders = re.findall(r"\{context\.(\w+)\}", path_template)
        
        resolved_path = path_template
        for key in placeholders:
            value = self.get(key)
            if value is None:
                raise ValueError(f"Could not resolve placeholder '{key}' in path template: '{path_template}'")
            
            # Replace placeholder with its value from context
            resolved_path = resolved_path.replace(f"{{context.{key}}}", str(value))
            
        return resolved_path
