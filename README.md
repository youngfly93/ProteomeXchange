# HLA Dataset Annotation Project

Automated annotation system for HLA immunopeptidomics datasets from proteomics repositories.

## Overview

This project provides automated classification of proteomics datasets by:
- **HLA Class**: I, II, I/II, or Unspecified
- **Analysis Scenario**: Cancer, Autoimmune, Infection, Normal, Mixed, Immunology
- **Disease Type**: Specific disease classifications (84 different types)

## Data

The repository contains metadata for **83 HLA datasets** from various sources:
- **79 PXD** (PRIDE) datasets
- **3 MSV** (MassIVE) datasets  
- **1 PASS** dataset

### Disease Distribution
- **Cancer studies**: 49 datasets (59%)
- **Autoimmune diseases**: 10 datasets (12%)
- **Normal/reference**: 8 datasets (10%)
- **Infections**: 6 datasets (7%)
- **Immunology**: 6 datasets (7%)
- **Mixed conditions**: 2 datasets (2%)

### HLA Class Distribution
- **Class I**: 58 datasets (70%)
- **Class II**: 11 datasets (13%)
- **Both I/II**: 11 datasets (13%)
- **Unspecified**: 3 datasets (4%)

## Files

| File | Description |
|------|-------------|
| `meta.txt` | Original metadata (Chinese/English mixed) |
| `dataset_annotation.tsv` | Standardized English annotations |
| `needs_manual.csv` | 8 cases requiring manual review |

## Configuration

- `diseases.yml`: Disease classification patterns
- `scenarios.yml`: Analysis scenario patterns  
- `hla_patterns.yml`: HLA class detection patterns

## Scripts

- `scripts/annotate.py`: Full API-based annotation (for new datasets)
- `scripts/process_existing.py`: Process existing meta.txt data
- `tests/test_classification.py`: Validation tests

## Usage

### Process Existing Data
```bash
python scripts/process_existing.py
```

### Run Tests
```bash
python -m pytest tests/test_classification.py -v
```

### Install Dependencies
```bash
pip install -r requirements.txt
# or
conda env create -f environment.yml
```

## Automation

GitHub Actions automatically:
1. Runs annotation when files change
2. Validates results with tests
3. Commits updated annotations

## Manual Review Cases

8 datasets require manual review (in `needs_manual.csv`):
- 3 with unspecified HLA class and scenario
- 5 normal/immunology samples with unspecified disease types

## Quality Validation

All annotations validated against known test cases:
- **PXD014397**: HLA Class I, Cancer, Melanoma ✓
- **PXD022930**: HLA Class II, Normal, Cell Line ✓  
- **PXD025499**: HLA I/II, Infection, COVID-19 ✓