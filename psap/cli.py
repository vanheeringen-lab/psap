"""Console script for psap."""
import argparse
import sys
from pathlib import Path
from .util import export_matrix
import psap


def main():
    """Console script for psap."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=psap.__version__)
    parser.add_argument(
        "-dbf",
        "--db_fasta",
        default=None,
        required=True,
        help='Path to proteome fasta file'
    )
    parser.add_argument(
        "-o",
        "--out",
        default="~",
        required=False,
        help='Output directory for serialized training-set'
    )
    parser.add_argument(
        "-cc",
        "--class_column",
        default=None,
        required=True,
        help="class column name for training-set"
    )
    args = parser.parse_args()
    basename = Path(args.db_fasta).stem
    # Pickle training-set
    export_matrix(basename, args.db_fasta, args.out, args.class_cloumn)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
