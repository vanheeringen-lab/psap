import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import datetime
from random import sample
from tqdm.notebook import tqdm
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import RandomForestClassifier
import sklearn_json as skljson
from pathlib import Path
from .matrix import MakeMatrix
from loguru import logger


def annotate(df, labels=None):
    if labels is None:
        labels = Path(__file__).parent / "data/assets/uniprot_ids.txt"
    with open(labels) as f:
        uniprot_ids = [line.rstrip() for line in f]
    logger.debug("Adding known llps class labels from {l}", l=labels)
    df["llps"] = 0
    for prot_id in uniprot_ids:
        df.loc[df["uniprot_id"] == prot_id, "llps"] = 1
        if (len(df.loc[df["uniprot_id"] == prot_id])) == 0:
            logger.info(prot_id + " is not found.")
    return df


def export_matrix(prefix="", fasta_path="", out_path=""):
    # Change pathing
    """Generates and saves a file which contains features of a protein sequence.
    Parameters:
        name: Name of the file.
        fasta_path: Path of the fasta file which needs to be featured.
    """
    logger.debug("Starting to add biochemical features to {f}", f=fasta_path)
    data = MakeMatrix(fasta_path)
    now = datetime.datetime.now()
    date = str(now.day) + "-" + str(now.month) + "-" + str(now.year)
    # Write data frame to csv
    out_dir = Path(out_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    data.df.to_csv(out_dir / f"{prefix}_psap_matrix_{date}.csv")
    return data


def preprocess_and_scaledata(data):
    """
    Wrapper for preprocess_data. Performs Min/Max scaling and centering of a dataset.
    ----------
    data : pandas DataFrame
            DataFrame with train/test instances
    Returns
    -------
    processed_data: Pandas DataFrame
                    Scaled and centered data
    """
    try:
        data = data.drop(["uniprot_id", "PRDaa", "HydroPhobicIndex"], axis=1)
    except KeyError:
        data = data.drop(["uniprot_id", "HydroPhobicIndex"], axis=1)
    logger.debug("Pre-processing and scaling dataset")
    data = data.fillna(value=0)
    scaler = MinMaxScaler()
    df = data.copy()
    processed_data = df.fillna(0)
    processed_data = preprocess_data(processed_data, scaler)
    # processed_data = remove_correlating_features(processed_data, cutoff=.95)
    # processed_data = remove_low_variance_features(processed_data, variance_cutoff=0.08)
    return processed_data


def preprocess_data(df, scaler):
    info = df.select_dtypes(include=["object"])
    X = df._get_numeric_data()
    columns = X.columns
    X = scaler.fit_transform(X)
    X = pd.DataFrame(X, columns=columns)
    X = X.merge(info, how="outer", left_index=True, right_index=True)
    return X


def train(
    path,
    prefix="",
    labels=None,
    out_dir="",
):
    """
    ----
    path: str
        Path to serialized/pickeld training set
    prefix:
        prefix for .csv file with prediction results and serialized model
    out_dir:
        path to create output folder.
    """
    mat = export_matrix(prefix=prefix, fasta_path=path, out_path=out_dir)
    data = annotate(mat.df, labels=labels)
    y = data["llps"]
    data_ps = preprocess_and_scaledata(data)
    # re-add class column after scaling
    data_ps["llps"] = y
    data_numeric = data_ps.select_dtypes([np.number])
    X = data_numeric.drop("llps", axis=1)
    y = data_numeric["llps"]
    # train random forest classifier
    logger.debug(
        "Training RF with {nin} instances and {nf} features",
        nf=len(X.columns),
        nin=len(X.index),
    )
    clf = RandomForestClassifier(
        n_jobs=32,
        class_weight="balanced",
        n_estimators=1200,
        criterion="entropy",
        random_state=42,
    )
    clf.fit(X, y)
    # write model to json
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"psap_model_{prefix}.json"
    logger.info("Writing trained RF classifier to {json}", json=out_file)
    skljson.to_json(clf, out_dir / f"psap_model_{prefix}.json")


def predict(
    path="",
    model=None,
    prefix="",
    out_dir="",
):
    """
    ----
    path: str
        Path to serialized/pickeld training set
    model: sklearn model
        Path to serialized RandomForest classifier (trained)
    prefix:
        prefix for .csv file with prediction results and serialized model
    out_dir:
        path to create output folder.
    """
    if model is None:
        model = Path(__file__).parent / "data/model/UP000005640_9606_llps.json"
    try:
        logger.info("Loading model: {m}", m=model)
        clf = skljson.from_json(model)
    except Exception:
        logger.error("classifier {mod} not found. Does the file exist?", mod=model)
    mat = export_matrix(prefix=prefix, fasta_path=path, out_path=out_dir)
    # Preprocessing
    data_ps = preprocess_and_scaledata(mat.df)
    X = data_ps.select_dtypes([np.number])
    logger.info("Predicting PSAP_score")
    psap_prediction = pd.DataFrame(index=data_ps["protein_name"])
    psap_prediction["PSAP_score"] = clf.predict_proba(X)[:, 1]
    psap_prediction["rank"] = 0
    rank = psap_prediction["PSAP_score"].rank(ascending=False)
    psap_prediction["rank"] = rank
    # # Make directory for output
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"prediction_{prefix}.csv"
    logger.info("Writing results to: {csv}", csv=out_file)
    psap_prediction.to_csv(out_file)
