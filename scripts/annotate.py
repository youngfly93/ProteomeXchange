#!/usr/bin/env python3
"""
HLA Dataset Annotation Script
Automatically classifies proteomics datasets by HLA class, analysis scenario, and disease type.
"""

import re
import csv
import yaml
import pandas as pd
import requests
import joblib
from pathlib import Path
from tqdm import tqdm
import time
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Setup cache
mem = joblib.Memory(".cache", verbose=0)

@mem.cache
def fetch_json(acc):
    """Fetch project and sample metadata from PRIDE API with caching."""
    base = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects"
    try:
        logger.info(f"Fetching data for {acc}")
        proj_response = requests.get(f"{base}/{acc}", timeout=15)
        proj_response.raise_for_status()
        proj = proj_response.json()
        
        samp_response = requests.get(f"{base}/{acc}/samples", timeout=15)
        samp_response.raise_for_status()
        samp = samp_response.json()
        
        # Add small delay to be respectful to API
        time.sleep(0.1)
        return proj, samp
    except requests.exceptions.RequestException as e:
        logger.warning(f"Failed to fetch data for {acc}: {e}")
        return None, None

def build_text(proj, samp):
    """Build searchable text from project and sample metadata."""
    if not proj:
        return ""
    
    fields = [
        proj.get("projectTitle", ""),
        proj.get("projectDescription", ""),
        " ".join(proj.get("keywords", [])),
        " ".join(attr.get("value", "") for attr in proj.get("additionalAttributes", [])),
    ]
    
    if samp:
        if isinstance(samp, list):
            for sample in samp:
                if isinstance(sample, dict):
                    # Handle different sample attribute formats
                    if "attributes" in sample:
                        if isinstance(sample["attributes"], str):
                            fields.append(sample["attributes"])
                        elif isinstance(sample["attributes"], list):
                            fields.extend([attr.get("value", "") for attr in sample["attributes"]])
                    
                    # Handle additional sample fields
                    for key in ["sampleDescription", "characteristics"]:
                        if key in sample:
                            if isinstance(sample[key], str):
                                fields.append(sample[key])
                            elif isinstance(sample[key], list):
                                fields.extend([str(item) for item in sample[key]])
    
    return " ".join(str(field) for field in fields if field).lower()

def load_patterns(file_path):
    """Load regex patterns from YAML file."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        logger.error(f"Pattern file not found: {file_path}")
        return {}

def classify(text, hla_patterns, scenario_patterns, disease_patterns):
    """Classify text based on regex patterns."""
    if not text:
        return "Unspecified", "Unspecified", "Unspecified"
    
    # 1) HLA Classification
    has_I = bool(re.search(hla_patterns.get("HLA_I", ""), text, re.IGNORECASE))
    has_II = bool(re.search(hla_patterns.get("HLA_II", ""), text, re.IGNORECASE))
    
    if has_I and has_II:
        hla = "I/II"
    elif has_I:
        hla = "I"
    elif has_II:
        hla = "II"
    else:
        hla = "Unspecified"
    
    # 2) Scenario Classification
    scenario_matches = []
    for scenario, pattern in scenario_patterns.items():
        if re.search(pattern, text, re.IGNORECASE):
            scenario_matches.append(scenario)
    
    if len(scenario_matches) > 1:
        scenario = "Mixed"
    elif len(scenario_matches) == 1:
        scenario = scenario_matches[0]
    else:
        scenario = "Unspecified"
    
    # 3) Disease Classification
    for disease, pattern in disease_patterns.items():
        if re.search(pattern, text, re.IGNORECASE):
            return hla, scenario, disease.replace("_", " ")
    
    return hla, scenario, "Unspecified"

def extract_accessions_from_meta(meta_file):
    """Extract accession numbers from existing meta.txt file."""
    accessions = []
    try:
        df = pd.read_csv(meta_file, sep='\t')
        accessions = df['all_accession'].tolist()
        logger.info(f"Extracted {len(accessions)} accessions from {meta_file}")
    except Exception as e:
        logger.error(f"Failed to read {meta_file}: {e}")
    
    return accessions

def main():
    """Main annotation function."""
    # Load configuration files
    hla_patterns = load_patterns("hla_patterns.yml")
    scenario_patterns = load_patterns("scenarios.yml")
    disease_patterns = load_patterns("diseases.yml")
    
    if not all([hla_patterns, scenario_patterns, disease_patterns]):
        logger.error("Failed to load one or more pattern files")
        return
    
    # Extract accessions from existing meta.txt
    accessions = extract_accessions_from_meta("meta.txt")
    
    if not accessions:
        logger.error("No accessions found to process")
        return
    
    # Results storage
    results = []
    needs_manual = []
    
    # Process each accession
    logger.info(f"Processing {len(accessions)} accessions...")
    
    for acc in tqdm(accessions, desc="Annotating datasets"):
        # Skip non-PXD accessions for now (PRIDE API specific)
        if not acc.startswith('PXD'):
            logger.info(f"Skipping non-PXD accession: {acc}")
            # Use existing data from meta.txt as fallback
            existing_df = pd.read_csv("meta.txt", sep='\t')
            existing_row = existing_df[existing_df['all_accession'] == acc]
            if not existing_row.empty:
                results.append({
                    'all_accession': acc,
                    'HLA': existing_row.iloc[0]['HLA(I/II)'],
                    'Scenario': existing_row.iloc[0]['分析场景'],
                    'Disease': existing_row.iloc[0]['疾病类型']
                })
            continue
        
        # Fetch metadata
        proj, samp = fetch_json(acc)
        text = build_text(proj, samp)
        
        # Classify
        hla, scenario, disease = classify(text, hla_patterns, scenario_patterns, disease_patterns)
        
        result = {
            'all_accession': acc,
            'HLA': hla,
            'Scenario': scenario,
            'Disease': disease
        }
        
        results.append(result)
        
        # Flag for manual review if any field is unspecified
        if hla == "Unspecified" or scenario == "Unspecified" or disease == "Unspecified":
            needs_manual.append(result)
    
    # Write results
    df_results = pd.DataFrame(results)
    df_results.to_csv("dataset_annotation.tsv", sep='\t', index=False)
    logger.info(f"Saved {len(results)} annotations to dataset_annotation.tsv")
    
    # Write manual review cases
    if needs_manual:
        df_manual = pd.DataFrame(needs_manual)
        df_manual.to_csv("needs_manual.csv", index=False)
        logger.info(f"Found {len(needs_manual)} cases requiring manual review in needs_manual.csv")
    
    # Summary statistics
    logger.info("=== Annotation Summary ===")
    logger.info(f"Total processed: {len(results)}")
    logger.info(f"HLA distribution:")
    for hla_type in df_results['HLA'].value_counts().items():
        logger.info(f"  {hla_type[0]}: {hla_type[1]}")
    logger.info(f"Scenario distribution:")
    for scenario in df_results['Scenario'].value_counts().items():
        logger.info(f"  {scenario[0]}: {scenario[1]}")
    logger.info(f"Manual review needed: {len(needs_manual)}")

if __name__ == "__main__":
    main()