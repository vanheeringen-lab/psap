"""Console script for psap."""
import argparse
import sys
from pathlib import Path
from psap.util import export_matrix
from psap.classifier import train_model, psap_predict


def main():
    """Console script for psap."""
    parser = argparse.ArgumentParser()
    # parser.add_argument("-v", "--version", action="version", version=psap.__version__)
    subparsers = parser.add_subparsers(dest="command")
    annotate = subparsers.add_parser(
        "annotate",
        help="adds biochemical features to a set of protein sequences in fasta format and writes it to a serialized data frame",
    )
    train = subparsers.add_parser("train", help="train psap model")
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
    train.add_argument(
        "-df",
        "--data_frame",
        default=None,
        required=True,
        help="annotated training set",
    )
    train.add_argument(
        "-out",
        "--out_dir",
        default=None,
        required=True,
        help="output directory for trained model",
    )
    pp.add_argument(
        "-df",
        "--data_frame",
        default=None,
        required=True,
        help="data frame with training data",
    )
    pp.add_argument(
        "-m",
        "--model",
        default=None,
        required=True,
        help="trained psap model",
    )
    pp.add_argument(
        "-out",
        "--out_dir",
        default=None,
        required=True,
        help=" serialized data frame with training data",
    )
    args = parser.parse_args()
    # Pickle training-set
    if args.command == "annotate":
        export_matrix(
            name=Path(args.db_fasta).stem, fasta_path=args.db_fasta, out_path=args.out
        )
    if args.command == "train":
        train_model(
            training_data=args.data_frame,
            prefix=Path(args.out_dir).stem,
            out_dir=args.out_dir,
        )
    elif args.command == "pp":
        psap_predict(
            test_data=args.data_frame,
            model=args.model,
            prefix=Path(args.out_dir).stem,
            out_dir=args.out_dir,
        )
    else:
        print("Incorrect subparser selected")


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
