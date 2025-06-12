#!/usr/bin/env python3
"""
HLA Dataset Annotation Script (Fixed Version)
Automatically classifies proteomics datasets by HLA class, analysis scenario, and disease type.
Focuses on project-level metadata to avoid API issues with sample endpoints.
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
def fetch_project_only(acc):
    """Fetch only project metadata from PRIDE API with caching."""
    base = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects"
    try:
        logger.info(f"Fetching project data for {acc}")
        proj_response = requests.get(f"{base}/{acc}", timeout=15)
        proj_response.raise_for_status()
        proj = proj_response.json()
        
        # Add small delay to be respectful to API
        time.sleep(0.1)
        return proj
    except requests.exceptions.RequestException as e:
        logger.warning(f"Failed to fetch project data for {acc}: {e}")
        return None

def build_text_from_project(proj):
    """Build searchable text from project metadata only."""
    if not proj:
        return ""
    
    def safe_join(items):
        """Safely join list items, handling dicts and other types."""
        if not items:
            return ""
        result = []
        for item in items:
            if isinstance(item, dict):
                # Extract text from dict values
                result.extend(str(v) for v in item.values() if v)
            elif isinstance(item, str):
                result.append(item)
            elif item:
                result.append(str(item))
        return " ".join(result)
    
    fields = [
        proj.get("projectTitle", ""),
        proj.get("projectDescription", ""),
        safe_join(proj.get("keywords", [])),
        " ".join(attr.get("value", "") for attr in proj.get("additionalAttributes", [])),
        safe_join(proj.get("instruments", [])),
        safe_join(proj.get("species", [])),
        safe_join(proj.get("tissues", [])),
        safe_join(proj.get("ptmList", [])),
        proj.get("doi", ""),
        str(proj.get("publicationDate", "")),
    ]
    
    # Also include submission details if available
    if "submissionType" in proj:
        fields.append(str(proj["submissionType"]))
    
    if "projectTags" in proj:
        fields.append(safe_join(proj.get("projectTags", [])))
    
    return " ".join(str(field) for field in fields if field).lower()

def load_patterns(file_path):
    """Load regex patterns from YAML file."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        logger.error(f"Pattern file not found: {file_path}")
        return {}

def enhanced_classify(text, hla_patterns, scenario_patterns, disease_patterns):
    """Enhanced classification with better pattern matching."""
    if not text:
        return "Unspecified", "Unspecified", "Unspecified"
    
    # Debug: print first few characters to see what we're working with
    logger.debug(f"Classification text sample: {text[:200]}...")
    
    # 1) HLA Classification - Enhanced patterns
    hla_i_patterns = [
        r'\bhla[- ]?i\b',
        r'\bclass[- ]?i\b',
        r'\bmhc[- ]?i\b',
        r'\bhla[- ]?class[- ]?i\b',
        r'\bmhc[- ]?class[- ]?i\b',
        r'\bhla-a\b', r'\bhla-b\b', r'\bhla-c\b',  # Specific HLA-I alleles
        r'\bh-2[dk]\b',  # Mouse HLA-I
    ]
    
    hla_ii_patterns = [
        r'\bhla[- ]?ii\b',
        r'\bclass[- ]?ii\b', 
        r'\bmhc[- ]?ii\b',
        r'\bhla[- ]?class[- ]?ii\b',
        r'\bmhc[- ]?class[- ]?ii\b',
        r'\bhla-dr\b', r'\bhla-dq\b', r'\bhla-dp\b',  # Specific HLA-II alleles
        r'\bh-2[ai]\b',  # Mouse HLA-II
    ]
    
    has_I = any(re.search(pattern, text, re.IGNORECASE) for pattern in hla_i_patterns)
    has_II = any(re.search(pattern, text, re.IGNORECASE) for pattern in hla_ii_patterns)
    
    if has_I and has_II:
        hla = "I/II"
    elif has_I:
        hla = "I"
    elif has_II:
        hla = "II"
    else:
        hla = "Unspecified"
    
    # 2) Scenario Classification - Enhanced
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
    
    # 3) Disease Classification - Enhanced with priority order
    # Check most specific diseases first
    disease_priority = [
        'COVID-19', 'Melanoma', 'Breast_Cancer', 'Lung_Cancer', 
        'Ovarian_Cancer', 'Hepatocellular_Carcinoma', 'Colorectal_Cancer',
        'Glioblastoma', 'Type_1_Diabetes', 'Multiple_Sclerosis',
        'Rheumatoid_Arthritis', 'Celiac_Disease', 'Behcets_Disease',
        'HIV_Infection', 'Tuberculosis', 'Influenza',
        'Acute_Myeloid_Leukemia', 'B_cell_Lymphoma', 'Chronic_Myeloid_Leukemia',
        'Cell_Line_Reference', 'Cancer', 'Autoimmune'
    ]
    
    for disease in disease_priority:
        if disease in disease_patterns:
            pattern = disease_patterns[disease]
            if re.search(pattern, text, re.IGNORECASE):
                return hla, scenario, disease.replace("_", " ")
    
    # Fallback to all patterns
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

def handle_non_pxd_accessions(acc, meta_file):
    """Handle non-PXD accessions using existing meta.txt data."""
    try:
        existing_df = pd.read_csv(meta_file, sep='\t')
        existing_row = existing_df[existing_df['all_accession'] == acc]
        if not existing_row.empty:
            return {
                'all_accession': acc,
                'HLA': existing_row.iloc[0]['HLA(I/II)'],
                'Scenario': existing_row.iloc[0]['分析场景'],
                'Disease': existing_row.iloc[0]['疾病类型']
            }
    except Exception as e:
        logger.error(f"Failed to get fallback data for {acc}: {e}")
    
    return {
        'all_accession': acc,
        'HLA': 'Unspecified',
        'Scenario': 'Unspecified', 
        'Disease': 'Unspecified'
    }

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
    api_success = 0
    api_failures = 0
    
    # Process each accession
    logger.info(f"Processing {len(accessions)} accessions...")
    
    for acc in tqdm(accessions, desc="Annotating datasets"):
        # Handle non-PXD accessions with existing data
        if not acc.startswith('PXD'):
            logger.info(f"Using fallback data for non-PXD accession: {acc}")
            result = handle_non_pxd_accessions(acc, "meta.txt")
            results.append(result)
            continue
        
        # Fetch project metadata only (no samples)
        proj = fetch_project_only(acc)
        
        if proj is None:
            logger.warning(f"Could not fetch project data for {acc}, using fallback")
            result = handle_non_pxd_accessions(acc, "meta.txt")
            api_failures += 1
        else:
            # Build text from project metadata only
            text = build_text_from_project(proj)
            
            if text:
                logger.debug(f"Text for {acc}: {text[:100]}...")
            
            # Classify using enhanced function
            hla, scenario, disease = enhanced_classify(text, hla_patterns, scenario_patterns, disease_patterns)
            
            result = {
                'all_accession': acc,
                'HLA': hla,
                'Scenario': scenario,
                'Disease': disease
            }
            api_success += 1
        
        results.append(result)
        
        # Flag for manual review if any field is unspecified
        if result['HLA'] == "Unspecified" or result['Scenario'] == "Unspecified" or result['Disease'] == "Unspecified":
            needs_manual.append(result)
    
    # Write results
    df_results = pd.DataFrame(results)
    df_results.to_csv("dataset_annotation_fixed.tsv", sep='\t', index=False)
    logger.info(f"Saved {len(results)} annotations to dataset_annotation_fixed.tsv")
    
    # Write manual review cases
    if needs_manual:
        df_manual = pd.DataFrame(needs_manual)
        df_manual.to_csv("needs_manual_fixed.csv", index=False)
        logger.info(f"Found {len(needs_manual)} cases requiring manual review in needs_manual_fixed.csv")
    
    # Summary statistics
    logger.info("=== Fixed Annotation Summary ===")
    logger.info(f"Total processed: {len(results)}")
    logger.info(f"API successes: {api_success}")
    logger.info(f"API failures: {api_failures}")
    logger.info(f"HLA distribution:")
    for hla_type in df_results['HLA'].value_counts().items():
        logger.info(f"  {hla_type[0]}: {hla_type[1]}")
    logger.info(f"Scenario distribution:")
    for scenario in df_results['Scenario'].value_counts().items():
        logger.info(f"  {scenario[0]}: {scenario[1]}")
    logger.info(f"Manual review needed: {len(needs_manual)}")

if __name__ == "__main__":
    main()