#!/usr/bin/env python

"""Tests for `psap` package."""


import pytest
import sklearn_json as skljson
from pathlib import Path
from psap.matrix import MakeMatrix


def test_model_exists():
    """Verify that the default model is part of the installation."""
    model_internal = (
        Path(__file__).parent.parent / "psap/data/model/UP000005640_9606_llps.json"
    )
    assert Path(model_internal).is_file()


def test_wl_exists():
    """Verify that the default whitelist is part of the installation."""
    wl_internal = Path(__file__).parent.parent / "psap/data/assets/uniprot_ids.txt"
    assert Path(wl_internal).is_file()


def test_fasta_exists():
    fasta = Path(__file__).parent.parent / "psap/data/testing/testset.fasta"
    assert Path(fasta).is_file()


def test_make_matrix():
    fasta = Path(__file__).parent.parent / "psap/data/testing/testset.fasta"
    matrix = MakeMatrix(fasta)
    assert len(matrix.df.columns) == 99


def test_feature_equal():
    fasta = Path(__file__).parent.parent / "psap/data/testing/testset.fasta"
    model_internal = (
        Path(__file__).parent.parent / "psap/data/model/UP000005640_9606_llps.json"
    )
    matrix = MakeMatrix(fasta)
    clf = skljson.from_json(model_internal)
    assert clf.n_features_ == len(matrix.df.columns) - 4
