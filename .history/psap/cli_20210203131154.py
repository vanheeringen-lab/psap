"""Console script for psap."""
import argparse
import sys
from pathlib import Path
from psap import fasta2df

def main():
    """Console script for psap."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="0.1.4-alpha"
    )
    parser.add_argument(
    "-dbf",
    "--db_fasta",
    )
    parser.add_argument(
    "-o",
    "--out",
    )
    args = parser.parse_args()
    basename = Path(args.db_fasta).stem
    #Export pickeled data frame
    fasta2df(basename, args.db_fasta, args.out_path)
    
if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
