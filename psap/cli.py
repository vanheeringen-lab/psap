"""Console script for psap."""
import argparse
import sys
from pathlib import Path
import psap
from psap.classifier import train, predict, cval, export_matrix


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
        help="adds biochemical features to a set of protein sequences in fasta format and writes it to a serialized data frame",
    )
    psap_train = subparsers.add_parser("train", help="train psap model")
    psap_predict = subparsers.add_parser("predict", help="predict classes")
    psap_cval = subparsers.add_parser(
        "cval", help="evaluate model using cross validation"
    )
    psap_annotate.add_argument(
        "-f",
        "--fasta",
        default=None,
        required=True,
        help="Path to peptide fasta file",
    )
    psap_annotate.add_argument(
        "-o",
        "--out",
        default="~",
        required=False,
        help="Output directory for annotated and serialized (pkl) data frame",
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
        default=None,
        required=True,
        help="Output directory for trained and serialized RandomForest classifier",
    )
    psap_train.add_argument(
        "-l",
        "--labels",
        default=Path(__file__).parent / "data/assets/uniprot_ids.txt",
        required=False,
        help="class labels for training instances",
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
        default=Path(__file__).parent / "data/model/UP000005640_9606_llps.json",
        required=False,
        help="Path to serialized RandomForest model",
    )
    psap_predict.add_argument(
        "-o",
        "--out",
        default=None,
        required=True,
        help="Output directory for prediction results",
    )
    psap_predict.add_argument(
        "-l",
        "--labels",
        default=Path(__file__).parent / "data/assets/uniprot_ids.txt",
        required=False,
        help="class labels for training instances",
    )
    psap_cval.add_argument(
        "-f",
        "--fasta",
        default=None,
        required=True,
        help="Path to peptide fasta file",
    )
    psap_cval.add_argument(
        "-o",
        "--out",
        default=None,
        required=True,
        help="Output directory for prediction results",
    )
    args = parser.parse_args()
    # Pickle training-set
    if args.command not in ["train", "predict", "annotate", "cval"]:
        parser.print_help()

    if args.command == "annotate":
        export_matrix(
            name=Path(args.fasta).stem, fasta_path=args.fasta, out_path=args.out
        )
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
            labels=args.labels,
            out_dir=args.out,
        )
    elif args.command == "cval":
        cval(
            path=args.data_frame,
            prefix=Path(args.out).stem,
            out_dir=args.out,
        )


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
