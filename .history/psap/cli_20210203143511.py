"""Console script for psap."""
import argparse
import sys
from pathlib import Path
from psap import export_matrix

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
    export_matrix(basename, args.db_fasta, args.out, True, 'llps')
    
if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
