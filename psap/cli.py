"""Console script for psap."""
import argparse
import sys
from pathlib import Path
from psap.util import export_matrix
from psap.classifier import train_model, psap_predict, eval_model


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
    predict = subparsers.add_parser("predict", help="predict classes")
    cval = subparsers.add_parser("cval", help="evaluate model using cross validation")
    annotate.add_argument(
        "-f",
        "--fasta",
        default=None,
        required=True,
        help="Path to peptide fasta file",
    )
    annotate.add_argument(
        "-o",
        "--out",
        default="~",
        required=False,
        help="Output directory for annotated and serialized (pkl) data frame",
    )
    train.add_argument(
        "-df",
        "--data_frame",
        default=None,
        required=True,
        help="Path to annotated and serialized data frame (output from annotate command)",
    )
    train.add_argument(
        "-o",
        "--out",
        default=None,
        required=True,
        help="Output directory for trained and serialized RandomForest classifier",
    )
    predict.add_argument(
        "-df",
        "--data_frame",
        default=None,
        required=True,
        help="Path to annotated and serialized data frame (output from annotate)",
    )
    predict.add_argument(
        "-m",
        "--model",
        default="data/model/UP000005640_9606_llps.joblid",
        required=False,
        help="Path to serialized RandomForest model",
    )
    predict.add_argument(
        "-o",
        "--out",
        default=None,
        required=True,
        help="Output directory for prediction results",
    )
    cval.add_argument(
        "-df",
        "--data_frame",
        default=None,
        required=True,
        help="Path to annotated and serialized data frame (output from annotate command)",
    )
    cval.add_argument(
        "-o",
        "--out",
        default=None,
        required=True,
        help="Output directory for prediction results",
    )
    args = parser.parse_args()
    # Pickle training-set
    if args.command == "annotate":
        export_matrix(
            name=Path(args.fasta).stem, fasta_path=args.fasta, out_path=args.out
        )
    elif args.command == "train":
        train_model(
            path=args.data_frame,
            prefix=Path(args.out).stem,
            out_dir=args.out,
        )
    elif args.command == "predict":
        psap_predict(
            path=args.data_frame,
            model=args.model,
            prefix=Path(args.out).stem,
            out_dir=args.out,
        )
    elif args.command == "cval":
        eval_model(
            path=args.data_frame,
            prefix=Path(args.out).stem,
            out_dir=args.out,
        )
    else:
        print("Incorrect subparser selected")


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
