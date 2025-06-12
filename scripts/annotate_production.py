#!/usr/bin/env python3
"""
HLA Dataset Annotation Script (Production Version)
Automatically classifies proteomics datasets by HLA class, analysis scenario, and disease type.
Enhanced with all fixes for production use.
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
from concurrent.futures import ThreadPoolExecutor, as_completed

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Setup cache
mem = joblib.Memory(".cache", verbose=0)

@mem.cache
def fetch_project_data(acc):
    """Fetch project metadata with v3/v2 fallback."""
    # Try v3 first (newer submissions)
    base_v3 = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects"
    base_v2 = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects"
    
    for base in [base_v3, base_v2]:
        try:
            logger.debug(f"Fetching {acc} from {base}")
            response = requests.get(f"{base}/{acc}", timeout=15)
            if response.status_code == 200:
                proj = response.json()
                time.sleep(0.1)  # Rate limiting
                return proj
        except requests.exceptions.RequestException as e:
            logger.debug(f"Failed to fetch from {base}: {e}")
            continue
    
    logger.warning(f"Failed to fetch project data for {acc} from both v3 and v2")
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
    
    if "submissionType" in proj:
        fields.append(str(proj["submissionType"]))
    
    if "projectTags" in proj:
        fields.append(safe_join(proj.get("projectTags", [])))
    
    return " ".join(str(field) for field in fields if field).lower()

def load_patterns(file_path):
    """Load regex patterns from YAML file with fallback."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        logger.warning(f"Pattern file not found: {file_path}, using defaults")
        # Return default patterns if file missing
        if "hla" in str(file_path).lower():
            return {
                "HLA_I": r"\b(hla[- ]?i|class[- ]?i|mhc[- ]?i)\b",
                "HLA_II": r"\b(hla[- ]?ii|class[- ]?ii|mhc[- ]?ii)\b"
            }
        elif "scenario" in str(file_path).lower():
            return {
                "Cancer": "(cancer|tumor|carcinoma|melanoma|neoplasm)",
                "Infection": "(virus|covid|hiv|influenza|bacteria|pathogen)",
                "Autoimmune": "(autoimmune|diabetes|celiac|lupus|arthritis)",
                "Normal": "(cell line|healthy|reference|normal control)",
                "Immunology": "(immunology|immune|immunological)",
                "Mixed": "(cancer.*virus|virus.*cancer)"
            }
        elif "disease" in str(file_path).lower():
            return {
                "COVID-19": "(covid[- ]?19|sars[- ]cov[- ]2)",
                "Melanoma": "(melanoma|skin cancer)",
                "Breast_Cancer": "(breast cancer|mammary carcinoma)",
                "Lung_Cancer": "(lung cancer|pulmonary carcinoma|nsclc|sclc)",
                "Type_1_Diabetes": "(type ?1 diabetes|t1d)",
                "Cancer": "(cancer|tumor|carcinoma|neoplasm)",
                "Unspecified": ""
            }
        return {}

def enhanced_classify(text, hla_patterns, scenario_patterns, disease_patterns):
    """Enhanced classification with better pattern matching."""
    if not text:
        return "Unspecified", "Unspecified", "Unspecified"
    
    # 1) HLA Classification - Using hardcoded patterns (more reliable)
    hla_i_patterns = [
        r'\bhla[- ]?i\b',
        r'\bclass[- ]?i\b',
        r'\bmhc[- ]?i\b',
        r'\bhla[- ]?class[- ]?i\b',
        r'\bmhc[- ]?class[- ]?i\b',
        r'\bhla-a\b', r'\bhla-b\b', r'\bhla-c\b',
        r'\bh-2[dk]\b',
    ]
    
    hla_ii_patterns = [
        r'\bhla[- ]?ii\b',
        r'\bclass[- ]?ii\b', 
        r'\bmhc[- ]?ii\b',
        r'\bhla[- ]?class[- ]?ii\b',
        r'\bmhc[- ]?class[- ]?ii\b',
        r'\bhla-dr\b', r'\bhla-dq\b', r'\bhla-dp\b',
        r'\bh-2[ai]\b',
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
    
    # 2) Scenario Classification
    scenario_matches = []
    for scenario, pattern in scenario_patterns.items():
        if pattern and re.search(pattern, text, re.IGNORECASE):
            scenario_matches.append(scenario)
    
    if len(scenario_matches) > 1:
        scenario = "Mixed"
    elif len(scenario_matches) == 1:
        scenario = scenario_matches[0]
    else:
        scenario = "Unspecified"
    
    # 3) Disease Classification with priority
    disease_priority = [
        'COVID-19', 'Melanoma', 'Breast_Cancer', 'Lung_Cancer', 
        'Type_1_Diabetes', 'Cancer', 'Unspecified'
    ]
    
    for disease in disease_priority:
        if disease in disease_patterns and disease_patterns[disease]:
            if re.search(disease_patterns[disease], text, re.IGNORECASE):
                return hla, scenario, disease.replace("_", " ")
    
    return hla, scenario, "Unspecified"

def extract_accessions_from_meta(meta_file):
    """Extract accession numbers from existing meta.txt file."""
    accessions = []
    try:
        df = pd.read_csv(meta_file, sep='\t')
        accessions = df.iloc[:, 0].tolist()  # First column regardless of name
        logger.info(f"Extracted {len(accessions)} accessions from {meta_file}")
    except Exception as e:
        logger.error(f"Failed to read {meta_file}: {e}")
    
    return accessions

def handle_non_pxd_accessions(acc, meta_file):
    """Handle non-PXD accessions with robust column access."""
    try:
        existing_df = pd.read_csv(meta_file, sep='\t')
        # Find row by first column (accession)
        existing_row = existing_df[existing_df.iloc[:, 0] == acc]
        
        if not existing_row.empty:
            row = existing_row.iloc[0]
            # Try multiple column name variations
            hla_names = ['HLA(I/II)', 'HLA', 'hla']
            scenario_names = ['分析场景', 'Scenario', 'scenario']
            disease_names = ['疾病类型', 'Disease', 'disease']
            
            hla = next((row.get(col, 'Unspecified') for col in hla_names if col in row), 'Unspecified')
            scenario = next((row.get(col, 'Unspecified') for col in scenario_names if col in row), 'Unspecified')
            disease = next((row.get(col, 'Unspecified') for col in disease_names if col in row), 'Unspecified')
            
            return {
                'all_accession': acc,
                'HLA(I/II)': hla,
                '分析场景': scenario,
                '疾病类型': disease
            }
    except Exception as e:
        logger.debug(f"Failed to get fallback data for {acc}: {e}")
    
    return {
        'all_accession': acc,
        'HLA(I/II)': 'Unspecified',
        '分析场景': 'Unspecified', 
        '疾病类型': 'Unspecified'
    }

def process_single_accession(acc, hla_patterns, scenario_patterns, disease_patterns, meta_file):
    """Process a single accession ID."""
    # Handle non-PXD accessions with existing data
    if not acc.startswith('PXD'):
        logger.info(f"Using fallback data for non-PXD accession: {acc}")
        return handle_non_pxd_accessions(acc, meta_file)
    
    # Fetch project metadata
    proj = fetch_project_data(acc)
    
    if proj is None:
        logger.warning(f"Could not fetch project data for {acc}, using fallback")
        return handle_non_pxd_accessions(acc, meta_file)
    
    # Build text and classify
    text = build_text_from_project(proj)
    hla, scenario, disease = enhanced_classify(text, hla_patterns, scenario_patterns, disease_patterns)
    
    return {
        'all_accession': acc,
        'HLA(I/II)': hla,
        '分析场景': scenario,
        '疾病类型': disease
    }

def main():
    """Main annotation function with parallel processing that preserves dataset order."""
    # Load configuration files
    hla_patterns = load_patterns("hla_patterns.yml")
    scenario_patterns = load_patterns("scenarios.yml")
    disease_patterns = load_patterns("diseases.yml")
    
    # Extract accessions
    meta_file = "meta.txt"
    accessions = extract_accessions_from_meta(meta_file)
    
    if not accessions:
        logger.error("No accessions found to process")
        return
    
    # Results storage - use dictionary to maintain order
    results_dict = {}
    needs_manual = []
    
    # Process with parallel execution
    logger.info(f"Processing {len(accessions)} accessions with parallel execution...")
    
    with ThreadPoolExecutor(max_workers=10) as executor:
        # Submit all tasks with index to preserve order
        future_to_info = {
            executor.submit(
                process_single_accession, 
                acc, hla_patterns, scenario_patterns, disease_patterns, meta_file
            ): (idx, acc) for idx, acc in enumerate(accessions)
        }
        
        # Process completed tasks with progress bar
        with tqdm(total=len(accessions), desc="Annotating datasets") as pbar:
            for future in as_completed(future_to_info):
                idx, acc = future_to_info[future]
                try:
                    result = future.result()
                    results_dict[idx] = result
                    
                    # Flag for manual review
                    if (result['HLA(I/II)'] == "Unspecified" or 
                        result['分析场景'] == "Unspecified" or 
                        result['疾病类型'] == "Unspecified"):
                        needs_manual.append(result)
                    
                except Exception as e:
                    logger.error(f"Error processing {acc}: {e}")
                    # Add failed result
                    results_dict[idx] = {
                        'all_accession': acc,
                        'HLA(I/II)': 'Unspecified',
                        '分析场景': 'Unspecified',
                        '疾病类型': 'Unspecified'
                    }
                
                pbar.update(1)
    
    # Convert to list preserving original order
    results = [results_dict[i] for i in sorted(results_dict.keys())]
    
    # Write results with correct column names
    df_results = pd.DataFrame(results)
    df_results.to_csv("dataset_annotation.tsv", sep='\t', index=False)
    logger.info(f"Saved {len(results)} annotations to dataset_annotation.tsv")
    
    # Write manual review cases
    if needs_manual:
        df_manual = pd.DataFrame(needs_manual)
        df_manual.to_csv("needs_manual.csv", index=False)
        logger.info(f"Found {len(needs_manual)} cases requiring manual review")
    
    # Summary statistics
    logger.info("=== Annotation Summary ===")
    logger.info(f"Total processed: {len(results)}")
    logger.info(f"HLA distribution:")
    for hla_type, count in df_results['HLA(I/II)'].value_counts().items():
        logger.info(f"  {hla_type}: {count}")
    logger.info(f"Scenario distribution:")
    for scenario, count in df_results['分析场景'].value_counts().items():
        logger.info(f"  {scenario}: {count}")
    logger.info(f"Disease distribution (top 10):")
    for disease, count in df_results['疾病类型'].value_counts().head(10).items():
        logger.info(f"  {disease}: {count}")
    logger.info(f"Manual review needed: {len(needs_manual)}")

if __name__ == "__main__":
    main()