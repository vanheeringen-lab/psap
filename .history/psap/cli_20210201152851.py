"""Console script for psap."""
import argparse
import sys


def main():
    """Console script for psap."""
    parser = argparse.ArgumentParser()
     )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="v1.0.1",
    )
    args = parser.parse_args()
    
if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
