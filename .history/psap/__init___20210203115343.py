"""Top-level package for psap."""
from .psap import MakeMatrix

__all__ = ["quantile_normalize"]

try:
    from .quantile_normalize import quantile_normalize_file  # noqa: F401

    __all__.append("quantile_normalize_file")
except ImportError:
    pass

__author__ = """Tilman Schaefers"""
__email__ = 'tilman.schaefers@ru.nl'
__version__ = '0.1.3'
