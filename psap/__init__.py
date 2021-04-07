"""Top-level package for psap."""
from .classifier import train, predict, export_matrix

__all__ = [
    "export_matrix",
    "train",
    "predict",
]


__author__ = """Tilman Schaefers"""
__email__ = "tilman.schaefers@ru.nl"

from ._version import get_versions

__version__ = get_versions()["version"]
##del get_versions
