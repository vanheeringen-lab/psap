"""Console script for psap."""
import argparse
import sys
from pathlib import Path
from psap.util import export_matrix
from psap.classifier import run_model


def main():
    """Console script for psap."""
    parser = argparse.ArgumentParser()
    # parser.add_argument("-v", "--version", action="version", version=psap.__version__)
    subparsers = parser.add_subparsers(dest="command")
    annotate = subparsers.add_parser(
        "annotate",
        help="adds biochemical features to a set of protein sequences in fasta format and writes it to a serialized data frame",
    )
    pp = subparsers.add_parser("pp", help="predict classes")
    annotate.add_argument(
        "-dbf",
        "--db_fasta",
        default=None,
        required=True,
        help="Path to proteome fasta file",
    )
    annotate.add_argument(
        "-o",
        "--out",
        default="~",
        required=False,
        help="Output directory for serialized training-set",
    )
    pp.add_argument(
        "-a",
        "--analysis",
        default=None,
        required=True,
        help="analysis name",
    )
    pp.add_argument(
        "-df",
        "--data_frame",
        default=None,
        required=True,
        help="data frame with training data",
    )
    pp.add_argument(
        "-out",
        "--out_dir",
        default=None,
        required=True,
        help=" serialized data frame with training data",
    )
    pp.add_argument(
        "-tdf",
        "--test_df",
        default=None,
        required=False,
        help="serialized data frame with test data",
    )
    args = parser.parse_args()
    # Pickle training-set
    if args.command == "annotate":
        export_matrix(Path(args.db_fasta).stem, args.db_fasta, args.out)
    elif args.command == "pp":
        run_model(args.analysis, args.data_frame, args.test_df, args.out_dir)
    else:
        print("Incorrect subparser selected")


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
