#!/usr/bin/env python3
"""
Unit tests for HLA dataset classification
"""

import pytest
import pandas as pd
import sys
import os

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from process_existing import standardize_hla, standardize_scenario, standardize_disease

class TestClassification:
    """Test classification functions."""
    
    def test_standardize_hla(self):
        """Test HLA standardization."""
        assert standardize_hla("I") == "I"
        assert standardize_hla("II") == "II"
        assert standardize_hla("I/II") == "I/II"
        assert standardize_hla("Unspecified") == "Unspecified"
        assert standardize_hla("invalid") == "Unspecified"
        assert standardize_hla("") == "Unspecified"
    
    def test_standardize_scenario(self):
        """Test scenario standardization."""
        assert standardize_scenario("Cancer") == "Cancer"
        assert standardize_scenario("Autoimmune") == "Autoimmune"
        assert standardize_scenario("Infection") == "Infection"
        assert standardize_scenario("Normal") == "Normal"
        assert standardize_scenario("Mixed") == "Mixed"
        assert standardize_scenario("Immunology") == "Immunology"
        assert standardize_scenario("invalid") == "Unspecified"
    
    def test_standardize_disease(self):
        """Test disease standardization."""
        assert standardize_disease("Cancer") == "Cancer"
        assert standardize_disease("Breast Cancer") == "Breast Cancer"
        assert standardize_disease("Behçet's Disease") == "Behcets Disease"
        assert standardize_disease("COVID-19") == "COVID-19"
        assert standardize_disease("Cell Line/Reference") == "Cell Line Reference"
        assert standardize_disease("Unspecified") == "Unspecified"

class TestOutputValidation:
    """Test output file validation against known cases."""
    
    @pytest.fixture
    def annotation_data(self):
        """Load annotation data for testing."""
        try:
            return pd.read_csv("dataset_annotation.tsv", sep='\t')
        except FileNotFoundError:
            pytest.skip("dataset_annotation.tsv not found")
    
    def test_known_cases(self, annotation_data):
        """Test specific known dataset classifications."""
        # Test PXD014397 - should be melanoma
        pxd014397 = annotation_data[annotation_data['all_accession'] == 'PXD014397']
        if not pxd014397.empty:
            row = pxd014397.iloc[0]
            assert row['HLA'] == 'I'
            assert row['Scenario'] == 'Cancer'
            assert row['Disease'] == 'Melanoma'
        
        # Test PXD022930 - should be HLA Class II, Normal, Cell Line
        pxd022930 = annotation_data[annotation_data['all_accession'] == 'PXD022930']
        if not pxd022930.empty:
            row = pxd022930.iloc[0]
            assert row['HLA'] == 'II'
            assert row['Scenario'] == 'Normal'
            assert row['Disease'] == 'Cell Line Reference'
        
        # Test PXD025499 - should be I/II, Infection, COVID-19
        pxd025499 = annotation_data[annotation_data['all_accession'] == 'PXD025499']
        if not pxd025499.empty:
            row = pxd025499.iloc[0]
            assert row['HLA'] == 'I/II'
            assert row['Scenario'] == 'Infection'
            assert row['Disease'] == 'COVID-19'
    
    def test_data_completeness(self, annotation_data):
        """Test that all records have required fields."""
        assert not annotation_data['all_accession'].isnull().any()
        assert not annotation_data['HLA'].isnull().any()
        assert not annotation_data['Scenario'].isnull().any()
        assert not annotation_data['Disease'].isnull().any()
    
    def test_hla_values(self, annotation_data):
        """Test that HLA values are within expected set."""
        valid_hla = {'I', 'II', 'I/II', 'Unspecified'}
        assert set(annotation_data['HLA'].unique()).issubset(valid_hla)
    
    def test_scenario_values(self, annotation_data):
        """Test that scenario values are within expected set."""
        valid_scenarios = {'Cancer', 'Autoimmune', 'Infection', 'Normal', 'Mixed', 'Immunology', 'Unspecified'}
        assert set(annotation_data['Scenario'].unique()).issubset(valid_scenarios)
    
    def test_accession_format(self, annotation_data):
        """Test that accession IDs follow expected format."""
        for acc in annotation_data['all_accession']:
            assert acc.startswith(('PXD', 'MSV', 'PASS')), f"Invalid accession format: {acc}"
    
    def test_record_count(self, annotation_data):
        """Test that we have the expected number of records."""
        # Should match the original meta.txt count
        assert len(annotation_data) == 83

if __name__ == "__main__":
    pytest.main([__file__, "-v"])