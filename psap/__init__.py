"""Top-level package for psap."""
from .util import export_matrix
from .psap import run_model

__all__ = ["export_matrix", "run_model"]


__author__ = """Tilman Schaefers"""
__email__ = "tilman.schaefers@ru.nl"

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
