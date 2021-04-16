#!/usr/bin/env python

"""Tests for `psap` package."""


import unittest
import sklearn_json as skljson
from pathlib import Path
from psap.matrix import MakeMatrix
class TestPsap(unittest.TestCase):
    """Tests for `psap` package."""

    def test_model_exists(self):
        """Verify that the default model is part of the installation."""
        model_internal = Path(__file__).parent.parent / "psap/data/model/UP000005640_9606_llps.json"
        self.assertTrue (Path(model_internal).is_file())
    
    def test_wl_exists(self):
        """Verify that the default whitelist is part of the installation."""
        wl_internal =  Path(__file__).parent.parent / "psap/data/assets/uniprot_ids.txt"
        self.assertTrue (Path(wl_internal).is_file())
    
    def test_fasta_exists(self):
        fasta =  Path(__file__).parent.parent / "psap/data/testing/testset.fasta"
        self.assertTrue (Path(fasta).is_file())
    
    def test_make_matrix(self):
        fasta =  Path(__file__).parent.parent / "psap/data/testing/testset.fasta"
        matrix = MakeMatrix(fasta)
        self.assertEqual (len(matrix.df.columns),91)
        
    def test_feature_equal(self):
        fasta =  Path(__file__).parent.parent / "psap/data/testing/testset.fasta"
        model_internal = Path(__file__).parent.parent / "psap/data/model/UP000005640_9606_llps.json"
        matrix = MakeMatrix(fasta)
        clf = skljson.from_json(model_internal)
        self.assertEqual (clf.n_features_,len(matrix.df.columns)-4, "different feature dimension between training/testset")   
       
    

    

            
        
        
        
        
        
