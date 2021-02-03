"""Console script for psap."""
import argparse
import sys


def main():
    """Console script for psap."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="0.1.3"
    )
    parser.add_argument(
    "-dbf",
    "--db_fasta",
    action="version",
    version="0.1.3"
    )

    
    args = parser.parse_args()
    
if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
