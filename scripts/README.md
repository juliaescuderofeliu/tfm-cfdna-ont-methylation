# Scripts overview and order of execution

This directory contains the scripts used for cfDNA methylation analysis in the Master’s Thesis.

Due to data governance restrictions, raw sequencing data (BAM, POD5, bedMethyl/CH3 files) are not included. The scripts document the analytical workflow and can be executed on appropriately formatted input data.

## Order of execution

### 1. Methylation calling (external / provider-generated)
Methylation calling was performed using `modkit pileup` on aligned BAM files.  
A representative shell script is provided for documentation purposes:

- `00_build_db/run_modkit_all.sh`

This step generates site-level methylation files (e.g. CH3 / bedMethyl).

---

### 2. Genome-wide window generation and aggregation
Scripts used to generate and process 1 Mb genomic windows:

- `00_build_db/windows1Mb_from_ch3.R`
- `00_build_db/merge_windows_1Mb.R`
- `00_build_db/label_and_qc_windows_1Mb.R`
- `00_build_db/filter_canonical_windows_1Mb.R`

These scripts aggregate CpG-level methylation into genome-wide 1 Mb windows and apply quality and canonical chromosome filtering.

---

### 3. Exploratory statistical analyses
Analyses performed on window-level and CpG-level data:

- `01_analysis/pca_windows_1Mb.R`
- `01_analysis/cpg_overlap_from_ch3.R`
- `01_analysis/cpg_matrix_pca_heatmap.R`
- `01_analysis/heatmap_topvar_windows.R`
- `01_analysis/diff_windows_1Mb.R`

These scripts generate PCA analyses, CpG overlap matrices, heatmaps, and exploratory differential methylation summaries.

---

### 4. Plot generation
Scripts used to generate figures included in the thesis:

- `02_plots/plot_delta_by_window.R`
- `02_plots/plot_diff_windows_1Mb.R`
- `02_plots/plot_top_delta_windows_1Mb_canonical.R`

---

### 5. Pilot analysis on chromosome 22
Scripts used for methodological validation on chromosome 22:

- `03_pilot_chr22/test70_chr22.R`

This pilot analysis was used to validate filtering thresholds and aggregation strategies before scaling to the whole genome.
