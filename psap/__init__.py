"""Top-level package for psap."""
from .util import export_matrix
from .classifier import train_model, psap_predict

__all__ = [
    "export_matrix",
    "train_model",
    "psap_predict",
    "eval_model",
]


__author__ = """Tilman Schaefers"""
__email__ = "tilman.schaefers@ru.nl"

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
