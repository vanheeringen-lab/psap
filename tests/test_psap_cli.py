#!/usr/bin/env python

"""Tests for `psap` package."""


import unittest
import sklearn_json as skljson
from pathlib import Path

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
    
    def test_model_features(self):
        """After importing the annotated data frame, the number of features should equal 87"""
        model_internal = Path(__file__).parent.parent / "psap/data/model/UP000005640_9606_llps.json"
        clf = skljson.from_json(model_internal)
        self.assertEqual (clf.n_features_,87)
    

            
        
        
        
        
        
