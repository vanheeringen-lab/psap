"""Console script for psap."""
import argparse
import sys
from pathlib import Path
import psap
from psap.classifier import annotate, train, predict, export_matrix
import os
from loguru import logger


def main():
    """Console script for psap."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"psap: v{psap.__version__}",
    )
    psap_annotate = subparsers.add_parser(
        "annotate",
        help="adds biochemical features to a set of protein sequences in fasta format and writes it to a csv file",
    )
    psap_train = subparsers.add_parser("train", help="train psap model")
    psap_predict = subparsers.add_parser("predict", help="predict PSAP score")
    psap_annotate.add_argument(
        "-f",
        "--fasta",
        default=None,
        required=True,
        help="Path to peptide fasta file",
    )
    psap_annotate.add_argument(
        "-l",
        "--labels",
        default=None,
        required=False,
        help=".txt file with llps uniprot ids (positive training labels)",
    )
    psap_annotate.add_argument(
        "-o",
        "--out",
        default=os.getcwd(),
        required=False,
        help="Output directory to store annotated data frame in .csv format",
    )
    psap_train.add_argument(
        "-f",
        "--fasta",
        default=None,
        required=True,
        help="Path to peptide fasta file",
    )
    psap_train.add_argument(
        "-o",
        "--out",
        default=os.getcwd(),
        required=False,
        help="Output directory to store trained RandomForest classifier in json format",
    )
    psap_train.add_argument(
        "-l",
        "--labels",
        default=None,
        required=False,
        help=".txt file with llps uniprot ids (positive training labels)",
    )
    psap_predict.add_argument(
        "-f",
        "--fasta",
        default=None,
        required=True,
        help="Path to peptide fasta file",
    )
    psap_predict.add_argument(
        "-m",
        "--model",
        required=False,
        help="Path to RandomForest model in json format",
    )
    psap_predict.add_argument(
        "-o",
        "--out",
        default=os.getcwd(),
        required=False,
        help="Output directory for psap prediction results",
    )
    psap_predict.add_argument(
        "-l",
        "--labels",
        default=None,
        required=False,
        help=".txt file with llps uniprot ids (positive training labels)",
    )
    args = parser.parse_args()
    # Pickle training-set
    if args.command not in ["train", "predict", "annotate"]:
        parser.print_help()

    if args.command == "annotate":
        mat = export_matrix(
            prefix=Path(args.fasta).stem,
            fasta_path=args.fasta,
            out_path=args.out,
        )
        annotate(mat.df, labels=args.labels)
    elif args.command == "train":
        train(
            path=args.fasta,
            prefix=Path(args.out).stem,
            labels=args.labels,
            out_dir=args.out,
        )
    elif args.command == "predict":
        predict(
            path=args.fasta,
            model=args.model,
            prefix=Path(args.out).stem,
            out_dir=args.out,
        )


if __name__ == "__main__":
    with logger.catch():
        sys.exit(main())  # pragma: no cover
