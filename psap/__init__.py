"""Top-level package for psap."""
from .util import export_matrix

__all__ = ["export_matrix"]


__author__ = """Tilman Schaefers"""
__email__ = "tilman.schaefers@ru.nl"
__version__ = "0.9.2a"

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
