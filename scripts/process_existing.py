#!/usr/bin/env python3
"""
Process existing meta.txt data to create standardized output format
"""

import pandas as pd
import yaml
import re
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_patterns(file_path):
    """Load regex patterns from YAML file."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        logger.error(f"Pattern file not found: {file_path}")
        return {}

def standardize_hla(hla_value):
    """Standardize HLA class notation."""
    hla_value = str(hla_value).strip()
    if hla_value in ["I", "II", "I/II", "Unspecified"]:
        return hla_value
    else:
        return "Unspecified"

def standardize_scenario(scenario_chinese):
    """Map Chinese scenario to English."""
    mapping = {
        "Cancer": "Cancer",
        "Autoimmune": "Autoimmune", 
        "Infection": "Infection",
        "Normal": "Normal",
        "Mixed": "Mixed",
        "Immunology": "Immunology",
        "Unspecified": "Unspecified"
    }
    return mapping.get(scenario_chinese, "Unspecified")

def standardize_disease(disease_chinese):
    """Map Chinese disease names to standardized English names."""
    disease_chinese = str(disease_chinese).strip()
    
    # Disease mapping dictionary
    disease_mapping = {
        "Cancer": "Cancer",
        "Breast Cancer": "Breast Cancer",
        "Ovarian Cancer": "Ovarian Cancer", 
        "Behçet's Disease": "Behcets Disease",
        "Melanoma": "Melanoma",
        "Cancer/Tumor": "Cancer",
        "Cell Line/Reference": "Cell Line Reference",
        "Hepatitis B Virus Infection, Hepatocellular Carcinoma": "Mixed Hepatitis B HCC",
        "COVID-19": "COVID-19",
        "Type 1 Diabetes": "Type 1 Diabetes",
        "Lung Cancer": "Lung Cancer",
        "Hepatocellular Carcinoma": "Hepatocellular Carcinoma",
        "Celiac Disease": "Celiac Disease",
        "Sarcoidosis": "Sarcoidosis",
        "Rheumatoid Arthritis, Lyme Arthritis": "Mixed Rheumatoid Lyme",
        "Glioblastoma": "Glioblastoma",
        "Immune-related Conditions": "Immune Related Conditions",
        "Mantle Cell Lymphoma": "Mantle Cell Lymphoma",
        "HIV/SIV Infection": "HIV Infection",
        "Multiple Sclerosis": "Multiple Sclerosis",
        "Birdshot Chorioretinopathy": "Birdshot Chorioretinopathy",
        "Ankylosing Spondylitis": "Ankylosing Spondylitis",
        "Colorectal Cancer": "Colorectal Cancer",
        "Lung Cancer, B-cell Acute Lymphoblastic Leukemia": "Mixed Lung Cancer B-ALL",
        "Meningioma": "Meningioma",
        "Chronic Myeloid Leukemia": "Chronic Myeloid Leukemia",
        "B-cell Lymphoma": "B-cell Lymphoma",
        "Acute Myeloid Leukemia": "Acute Myeloid Leukemia",
        "Tuberculosis": "Tuberculosis",
        "Influenza Virus Infection": "Influenza",
        "Diffuse Large B-cell Lymphoma": "Diffuse Large B-cell Lymphoma",
        "Melanoma, Lung Cancer": "Mixed Melanoma Lung Cancer",
        "Chronic Lymphocytic Leukemia": "Chronic Lymphocytic Leukemia",
        "Unspecified": "Unspecified"
    }
    
    return disease_mapping.get(disease_chinese, disease_chinese)

def main():
    """Process existing meta.txt data."""
    try:
        # Read existing data
        df = pd.read_csv("meta.txt", sep='\t')
        logger.info(f"Loaded {len(df)} records from meta.txt")
        
        # Standardize columns
        df_standardized = pd.DataFrame()
        df_standardized['all_accession'] = df['all_accession']
        df_standardized['HLA'] = df['HLA(I/II)'].apply(standardize_hla)
        df_standardized['Scenario'] = df['分析场景'].apply(standardize_scenario)
        df_standardized['Disease'] = df['疾病类型'].apply(standardize_disease)
        
        # Save standardized data
        df_standardized.to_csv("dataset_annotation.tsv", sep='\t', index=False)
        logger.info(f"Saved standardized annotations to dataset_annotation.tsv")
        
        # Identify cases needing manual review
        needs_manual = df_standardized[
            (df_standardized['HLA'] == 'Unspecified') |
            (df_standardized['Scenario'] == 'Unspecified') |
            (df_standardized['Disease'] == 'Unspecified')
        ]
        
        if len(needs_manual) > 0:
            needs_manual.to_csv("needs_manual.csv", index=False)
            logger.info(f"Found {len(needs_manual)} cases requiring manual review")
        
        # Generate summary statistics
        logger.info("=== Dataset Summary ===")
        logger.info(f"Total datasets: {len(df_standardized)}")
        
        logger.info("\nHLA Class distribution:")
        for hla_type, count in df_standardized['HLA'].value_counts().items():
            logger.info(f"  {hla_type}: {count}")
        
        logger.info("\nScenario distribution:")
        for scenario, count in df_standardized['Scenario'].value_counts().items():
            logger.info(f"  {scenario}: {count}")
        
        logger.info("\nTop 10 Disease types:")
        for disease, count in df_standardized['Disease'].value_counts().head(10).items():
            logger.info(f"  {disease}: {count}")
        
        logger.info(f"\nCases needing manual review: {len(needs_manual)}")
        
        # Dataset source distribution
        logger.info("\nDataset source distribution:")
        source_counts = df_standardized['all_accession'].str.extract(r'^([A-Z]+)')[0].value_counts()
        for source, count in source_counts.items():
            logger.info(f"  {source}: {count}")
            
    except Exception as e:
        logger.error(f"Error processing data: {e}")
        raise

if __name__ == "__main__":
    main()