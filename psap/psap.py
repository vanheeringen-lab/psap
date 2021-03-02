"""Main module."""
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from tqdm.auto import tqdm
import time
from scipy import signal
import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from random import sample
from tqdm.notebook import tqdm
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import RandomForestClassifier

RESIDUES = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]

# Kyte & Doolittle {kd} index of hydrophobicity
HP = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
    "U": 0.0,
}


class MakeMatrix:
    def __init__(self, dbfasta):
        self.df = pd.DataFrame()
        self.dbfasta = dbfasta
        self.add_features()

    def add_features(self):
        executables = [
            "self.fasta2df()",
            "self.amino_acid_analysis()",
            "self.add_biochemical_combinations()",
            "self.add_lowcomplexity_features()",
        ]
        for e in executables:
            start = time.time()
            print(e)
            exec(e)
            end = time.time()
            print(str(round(end - start, 2)) + "s " + e)

    def fasta2df(self):
        rows = list()
        with open(self.dbfasta) as f:
            for record in SeqIO.parse(self.dbfasta, "fasta"):
                seqdict = dict()
                seq = str(record.seq)
                id = record.description.split("|")
                if id[0] == "sp":
                    uniprot_id = id[1]
                    name = id[2].split(" ")[0]
                    rows.append([name, uniprot_id, seq])
                elif id[0] == "tr":
                    uniprot_id = id[1]
                    name = id[2].split(" ")[0]
                    rows.append([name, uniprot_id, seq])
                else:
                    uniprot_id = id[0]
                    name = id[2].split(" ")[0]
                    rows.append([name, uniprot_id, seq])
        self.df = pd.DataFrame(rows, columns=["protein_name", "uniprot_id", "sequence"])

    def hydrophobic(self):
        for index, row in self.df.iterrows():
            hpilst = pd.Series(list(row["sequence"])).map(HP).tolist()
            self.df.loc[index, "HydroPhobicIndex"] = HydroPhobicIndex(hpilst)

    def amino_acid_analysis(self):
        for res in RESIDUES:
            self.df["fraction_" + res] = (
                self.df["sequence"].str.count(res) / self.df["sequence"].str.len()
            )
        self.df["length"] = self.df["sequence"].str.len()
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0]):
            # for index, row in self.df.iterrows():
            seq = row["sequence"]
            seqanalysis = ProteinAnalysis(seq)
            acidist = seqanalysis.get_amino_acids_percent()
            self.df.loc[index, "IEP"] = seqanalysis.isoelectric_point()
            if "X" not in seq and "B" not in seq:
                self.df.loc[index, "molecular_weight"] = seqanalysis.molecular_weight()
            if "U" not in seq and "X" not in seq and "B" not in seq:
                self.df.loc[index, "gravy"] = seqanalysis.gravy()

    def add_iupred_features(self):
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0]):
            # for index, row in self.df.iterrows():
            idr = row["iupred"].glob[0]
            self.df.loc[index, "idr_percetage"] = sum(i > 0.5 for i in list(idr))
            self.df.loc[index, "idr_50"] = sum(i > 0.5 for i in list(idr)) / len(
                str(row["sequence"])
            )
            self.df.loc[index, "idr_60"] = sum(i > 0.6 for i in list(idr)) / len(
                str(row["sequence"])
            )
            self.df.loc[index, "idr_70"] = sum(i > 0.7 for i in list(idr)) / len(
                str(row["sequence"])
            )
            self.df.loc[index, "idr_80"] = sum(i > 0.8 for i in list(idr)) / len(
                str(row["sequence"])
            )
            self.df.loc[index, "idr_90"] = sum(i > 0.9 for i in list(idr)) / len(
                str(row["sequence"])
            )

    @staticmethod
    def convolve_signal(sig, window=25):
        win = signal.hann(window)
        sig = signal.convolve(sig, win, mode="same") / sum(win)
        return sig

    def add_hydrophobic_features(self):
        hpi0, hpi1, hpi2, hpi3, hpi4, hpi5 = (
            list(),
            list(),
            list(),
            list(),
            list(),
            list(),
        )
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0]):
            # for index, row in self.df.iterrows():
            sw = convolve_signal(row["HydroPhobicIndex"].hpilist, window=30)
            hpi0.append(sum(i < -1.5 for i in sw) / len(sw))
            # self.df.loc[index, 'hpi_<-1.5_frac'] = hpi
            hpi1.append(sum(i < -2.0 for i in sw) / len(sw))
            # self.df.loc[index, 'hpi_<-2.0_frac'] = hpi
            hpi2.append(sum(i < -2.5 for i in sw) / len(sw))
            # self.df.loc[index, 'hpi_<-2.5_frac'] = hpi
            hpi3.append(sum(i < -1.5 for i in sw))
            # self.df.loc[index, 'hpi_<-1.5'] = hpi
            hpi4.append(sum(i < -2.0 for i in sw))
            # self.df.loc[index, 'hpi_<-2.0'] = hpi
            hpi5.append(sum(i < -2.5 for i in sw))
            # self.df.loc[index, 'hpi_<-2.5'] = hpi
        self.df["hpi_<-1.5_frac"] = hpi0
        self.df["hpi_<-2.0_frac"] = hpi1
        self.df["hpi_<-2.5_frac"] = hpi2
        self.df["hpi_<-1.5"] = hpi3
        self.df["hpi_<-2.0"] = hpi4
        self.df["hpi_<-2.5"] = hpi5

    def add_biochemical_combinations(self):
        df = self.df
        df = df.assign(Asx=df["fraction_D"] + df["fraction_N"])
        df = df.assign(Glx=df["fraction_E"] + df["fraction_Q"])
        df = df.assign(Xle=df["fraction_I"] + df["fraction_L"])
        df = df.assign(
            Pos_charge=df["fraction_K"] + df["fraction_R"] + df["fraction_H"]
        )
        df = df.assign(Neg_charge=df["fraction_D"] + df["fraction_E"])
        df = df.assign(
            Aromatic=df["fraction_F"]
            + df["fraction_W"]
            + df["fraction_Y"]
            + df["fraction_H"]
        )
        df = df.assign(
            Alipatic=df["fraction_V"]
            + df["fraction_I"]
            + df["fraction_L"]
            + df["fraction_M"]
        )
        df = df.assign(
            Small=df["fraction_P"]
            + df["fraction_G"]
            + df["fraction_A"]
            + df["fraction_S"]
        )
        df = df.assign(
            Hydrophilic=(
                df["fraction_S"]
                + df["fraction_T"]
                + df["fraction_H"]
                + df["fraction_N"]
                + df["fraction_Q"]
                + df["fraction_E"]
                + df["fraction_D"]
                + df["fraction_K"]
                + df["fraction_R"]
            )
        )
        df = df.assign(
            Hydrophobic=(
                df["fraction_V"]
                + df["fraction_I"]
                + df["fraction_L"]
                + df["fraction_F"]
                + df["fraction_W"]
                + df["fraction_Y"]
                + df["fraction_M"]
            )
        )
        # Added in version 2
        for dimer in [
            "GV",
            "VG",
            "VP",
            "PG",
            "FG",
            "RG",
            "GR",
            "GG",
            "YG",
            "GS",
            "SG",
            "GA",
            "GF",
            "GD",
            "DS",
        ]:
            self.df[dimer] = self.df["sequence"].str.count(dimer)
        df = df.assign(
            alpha_helix=df["fraction_V"]
            + df["fraction_I"]
            + df["fraction_Y"]
            + df["fraction_F"]
            + df["fraction_W"]
            + df["fraction_L"]
        )
        df = df.assign(
            beta_turn=df["fraction_N"]
            + df["fraction_P"]
            + df["fraction_G"]
            + df["fraction_S"]
        )
        df = df.assign(
            beta_sheet=df["fraction_E"]
            + df["fraction_M"]
            + df["fraction_A"]
            + df["fraction_L"]
        )
        # Calculates the aromaticity value of a protein according to Lobry, 1994.
        # It is simply the relative frequency of Phe+Trp+Tyr.
        df = df.assign(
            aromaticity=df["fraction_F"] + df["fraction_W"] + df["fraction_Y"]
        )
        self.df = df
        del df

    def add_lowcomplexityscore(self):
        lcs_window = 20
        lcs_cutoff = 7
        for index, row in self.df.iterrows():
            seq = str(row["sequence"])
            if len(seq) > lcs_window + 1:
                sig = list()
                for i in range(len(seq)):
                    window = seq[i : i + lcs_window]
                    if len(window) == lcs_window:
                        acid_comp = len(list(set(window)))
                        sig.append(acid_comp)
                score = sum([1 if i <= 7 else 0 for i in sig])
                self.df.loc[index, "lcs_score"] = score
                self.df.loc[index, "lcs_fraction"] = score / len(sig)

    def add_lowcomplexity_features(self):
        n_window = 20
        cutoff = 7
        n_halfwindow = int(n_window / 2)
        lcs_lowest_complexity = list()
        lcs_scores = list()
        lcs_fractions = list()
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0]):
            # for index, row in self.df.iterrows():
            # Determine low complexity scores
            seq = str(row["sequence"])
            lcs_acids = list()
            sig = list()
            # New
            lc_bool = [False] * len(seq)
            for i in range(len(seq)):
                if i < n_halfwindow:
                    peptide = seq[:n_window]
                elif i + n_halfwindow > int(len(seq)):
                    peptide = seq[-n_window:]
                else:
                    peptide = seq[i - n_halfwindow : i + n_halfwindow]
                complexity = len(set(peptide))
                if complexity <= 7:
                    for bool_index in (i - n_halfwindow, i + n_halfwindow):
                        try:
                            lc_bool[bool_index] = True
                        except IndexError:
                            pass
                    lcs_acids.append(seq[i])
                sig.append(complexity)
            # Adding low complexity scores to list
            low_complexity_list = pd.DataFrame(
                {"bool": lc_bool, "acid": list(seq)}, index=None
            )
            lcs_lowest_complexity.append(min(sig))
            lcs_scores.append(
                len(low_complexity_list.loc[low_complexity_list["bool"] == True])
            )
            lcs_fractions.append(
                len(low_complexity_list.loc[low_complexity_list["bool"] == True])
                / len(seq)
            )
            low_complexity_list = pd.DataFrame(
                {"bool": lc_bool, "acid": list(seq)}, index=None
            )
            if len(lcs_acids) >= n_window:
                for i in RESIDUES:
                    self.df.loc[index, i + "_lcscore"] = len(
                        low_complexity_list.loc[
                            (low_complexity_list["bool"] == True)
                            & (low_complexity_list["acid"] == i)
                        ]
                    )
                    self.df.loc[index, i + "_lcfraction"] = len(
                        low_complexity_list.loc[
                            (low_complexity_list["bool"] == True)
                            & (low_complexity_list["acid"] == i)
                        ]
                    ) / len(lcs_acids)
        self.df["lcs_fractions"] = lcs_fractions
        self.df["lcs_scores"] = lcs_scores
        self.df["lcs_lowest_complexity"] = lcs_lowest_complexity


def preprocess_and_scaledata(data, instance):
    # try:
    #    data = data.drop(['iupred', 'HydroPhobicIndex', 'uniprot_id', 'PRDaa'], axis=1)
    # except KeyError:
    #    data = data.drop(['iupred', 'HydroPhobicIndex', 'uniprot_id'], axis=1)
    data = data.fillna(value=0)
    print(data.shape)
    print(
        "Number of phase separating proteins in dataset: "
        + str(data.loc[data[instance] == 1].shape[0])
    )
    print(data.shape)
    scaler = MinMaxScaler()
    df = data.copy()
    processed_data = df.fillna(0)
    print(processed_data)
    processed_data = preprocess_data(processed_data, scaler)
    # processed_data = remove_correlating_features(processed_data, cutoff=.95)
    # processed_data = remove_low_variance_features(processed_data, variance_cutoff=0.08)
    return processed_data


def preprocess_data(df, scaler):
    instance = "sample8"
    info = df.select_dtypes(include=["object"])
    y = df[instance]
    X = df.drop([instance], axis=1)
    X = X._get_numeric_data()
    columns = X.columns
    X = scaler.fit_transform(X)
    X = pd.DataFrame(X, columns=columns)
    X[instance] = y
    X = X.merge(info, how="outer", left_index=True, right_index=True)
    return X


def get_test_train_indexes(data, label, ratio=1, randomized=False):
    """
    Function: Will oversample the positive data with randomly selected negative samples.
    Returns: List with indexes which contain positive and negative samples.
    """
    positive_instances = set(data.loc[data[label] == 1].index)
    negative_instances = set(data.loc[data[label] == 0].index)
    n_positives = data.loc[data[label] == 1].shape[0]
    indexes = list()
    while len(negative_instances) >= 1:
        if len(negative_instances) > n_positives * ratio:
            sample_set = sample(negative_instances, (n_positives * ratio))
        else:
            sample_set = list(negative_instances)
            size = len(sample_set)
            short = (len(positive_instances) * ratio) - size
            shortage = sample(set(data.loc[data[label] == 0].index), short)
            sample_set = sample_set + shortage
        indexes.append((list(positive_instances) + list(sample_set)))
        negative_instances.difference_update(set(sample_set))
    return indexes


def predict_proteome(
    df,
    clf,
    instance,
    testing_size,
    feature_imp=True,
    remove_training=False,
    second_df=pd.DataFrame(),
):
    pd.set_option("mode.chained_assignment", None)
    prediction = df.select_dtypes(include="object")
    df = df.select_dtypes([np.number])
    if len(second_df) > 0:
        prediction = second_df.select_dtypes(include="object")
        second_df = second_df.select_dtypes([np.number])
    indexes = get_test_train_indexes(df, instance)
    count = 0
    fi_data = None
    for index in tqdm(indexes):
        df_fraction = df.loc[index]
        # Also consider X_test index for prediction in the proteome
        X = df_fraction.drop(instance, axis=1)
        y = df_fraction[instance]
        clf.fit(X, y)
        # Feature importance
        if feature_imp:
            # X = df_fraction.drop('llps', axis=1)
            fi = clf.feature_importances_
            fi = pd.DataFrame(fi).transpose()
            fi.columns = X.columns
            fi = fi.melt()
            fi = fi.sort_values("value", ascending=False)
            if fi_data is None:
                fi_data = fi
            else:
                fi_data = pd.merge(fi_data, fi, on="variable")
        # Make prediction
        if len(second_df) > 0:
            probability = clf.predict_proba(second_df.drop(instance, axis=1))[:, 1]
            prediction["probability_" + str(count)] = probability
        else:
            probability = clf.predict_proba(df.drop(instance, axis=1))[:, 1]
            prediction["probability_" + str(count)] = probability
        # Removing prediction that were used in the train test set.
        if remove_training:
            for i in index:
                prediction.loc[i, "probability_" + str(count)] = np.nan
        count += 1
    if feature_imp:
        return prediction, fi_data
    else:
        return prediction


def main(
    data,
    instance,
    ANALYSIS_NAME,
    second_df=pd.DataFrame(),
    second_df_bool=False,
    out_dir="/Users/tilman_work",
):
    # Make directory for output.
    try:
        os.mkdir(f"{out_dir}/{ANALYSIS_NAME}")
    except:
        print(
            f"Directory {ANALYSIS_NAME} already exists. Please choose another analysis name, or remove the directory {ANALYSIS_NAME}."
        )
    # Make prediction with random forest
    clf = RandomForestClassifier(max_depth=12, n_estimators=100)
    prediction, fi_data = predict_proteome(
        data, clf, instance=instance, testing_size=0.2, remove_training=False
    )
    # Predict second dataset
    # if second_df_bool:
    #    predict_other_df(clf, data, second_df, ANALYSIS_NAME)
    # Get Feature Importance
    # plot_feature_importance(fi_data, ANALYSIS_NAME)
    # Save prediction to .csv
    prediction.to_csv(f"{out_dir}/{ANALYSIS_NAME}/preidction_{ANALYSIS_NAME}.csv")


def run_model(ANALYSIS_NAME, instance, path):
    data = pd.read_pickle(path)
    data = preprocess_and_scaledata(data, instance)
    print("preprocessed and scaled dataset")
    main(data, instance, ANALYSIS_NAME)
