# cfDNA methylation analysis using Oxford Nanopore sequencing

This repository contains the bioinformatic scripts used in the Master’s Thesis:

**Oxford Nanopore cfDNA methylation in metastatic colorectal cancer and integration with precomputed CNAs**

## Overview

The pipeline processes cfDNA methylation data derived from Oxford Nanopore sequencing and performs:
- Quality filtering of methylation calls
- Genome-wide aggregation using 1 Mb genomic windows
- CpG-level overlap analysis across samples
- Exploratory statistical analyses (PCA, heatmaps)
- Qualitative integration with precomputed copy-number alteration (CNA) profiles

## Repository structure

- `scripts/`: R scripts used for methylation processing and analysis  
- `metadata/`: Sample annotation files  
- `results/`: Output tables and figures (not included due to data governance constraints)  
- `docs/`: Additional documentation

## Requirements

- R version 4.5.1
- MethylSeqR version 0.5.0
- Additional R packages: tidyverse, data.table, ggplot2, pheatmap

## Data availability

Raw sequencing data and methylation files are not publicly available due to patient privacy and data governance restrictions. The scripts provided here allow reproduction of the analytical workflow when applied to appropriately formatted input data.

## Author

Julia Escudero  
Master’s Thesis – Universitat Oberta de Catalunya (UOC)
