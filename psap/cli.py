"""Console script for psap."""
import argparse
import sys
from pathlib import Path
from .util import export_matrix
from .psap import run_model
import psap


def main():
    """Console script for psap."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", action="version", version=psap.__version__)
    subparsers = parser.add_subparsers(dest="command")
    train = subparsers.add_parser("train", help="create psap training-set")
    classify = subparsers.add_parser("classify", help="predict classes")
    train.add_argument(
        "-dbf",
        "--db_fasta",
        default=None,
        required=True,
        help="Path to proteome fasta file",
    )
    train.add_argument(
        "-o",
        "--out",
        default="~",
        required=False,
        help="Output directory for serialized training-set",
    )
    train.add_argument(
        "-cc",
        "--class_column",
        default=None,
        required=True,
        help="class column name for training-set",
    )
    classify.add_argument(
        "-a",
        "--analysis",
        default=None,
        required=True,
        help="analysis name",
    )
    classify.add_argument(
        "-i",
        "--instance",
        default=None,
        required=True,
        help="instance name",
    )
    classify.add_argument(
        "-df",
        "--data_frame",
        default=None,
        required=True,
        help="data frame with training data",
    )
    args = parser.parse_args()
    # Pickle training-set
    if args.command == "train":
        export_matrix(
            Path(args.db_fasta).stem, args.db_fasta, args.out, args.class_column
        )
    elif args.command == "classify":
        run_model(args.analysis, args.instance, args.data_frame)
    else:
        print("Incorrect subparser selected")


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
