"""Console script for psap."""
import argparse
import sys
import psap

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
    print(MakeMatrix(args.db_fasta))
    
if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
