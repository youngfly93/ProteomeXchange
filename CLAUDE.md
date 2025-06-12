# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a data repository containing HLA (Human Leukocyte Antigen) dataset metadata for immunopeptidomics research. The repository currently contains a single data file with proteomics dataset information.

## Data Structure

- `meta.txt`: Tab-separated metadata file containing proteomics dataset accessions with the following columns:
  - `all_accession`: Dataset accession numbers (PXD, MSV, PASS formats)
  - `HLA(I/II)`: HLA class type (I, II, I/II, or Unspecified)
  - `分析场景`: Analysis scenario (Cancer, Autoimmune, Infection, Normal, Mixed, Immunology)
  - `疾病类型`: Disease type (specific cancer types, autoimmune diseases, infections, etc.)

## Working with the Data

The metadata covers 84 datasets spanning various disease contexts including:
- Cancer studies (multiple cancer types)
- Autoimmune diseases
- Infectious diseases (COVID-19, HIV, Tuberculosis, etc.)
- Normal/reference samples
- Cell line studies

When analyzing or processing this data, consider the mixed language content (English headers with Chinese category descriptions) and the variety of dataset sources and disease classifications.