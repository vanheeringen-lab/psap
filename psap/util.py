from .matrix import MakeMatrix
import datetime
from pathlib import Path


def annotate(df, identifier_name):
    with open("data/assets/uniprot_ids.txt") as f:
        uniprot_ids = [line.rstrip() for line in f]
    print(uniprot_ids)
    df[identifier_name] = 0
    for prot_id in uniprot_ids:
        df.loc[df["uniprot_id"] == prot_id, identifier_name] = 1
        if (len(df.loc[df["uniprot_id"] == prot_id])) == 0:
            print(prot_id + " is not found.")
    return df


def export_matrix(name, fasta_path, out_path):
    # Change pathing
    """Generates and saves a file which contains features of a protein sequence.
    Parameters:
        name: Name of the file.
        fasta_path: Path of the fasta file which needs to be featured.
    """
    class_col = "llps"
    data = MakeMatrix(fasta_path)
    now = datetime.datetime.now()
    date = str(now.day) + "-" + str(now.month) + "-" + str(now.year)
    print("Adding labels to df")
    df_ann = annotate(data.df, class_col)
    # Write data frame to csv
    out_dir = Path(out_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    df_ann.to_csv(out_dir / f"{name}_{class_col}_{date}_ann.csv")
    return df_ann
