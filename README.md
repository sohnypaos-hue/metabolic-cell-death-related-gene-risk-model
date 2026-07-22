# A generalizable multi-pathway metabolic cell death signature predicts clinical outcomes and tumor microenvironmental states in bladder cancer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21491997.svg)](https://doi.org/10.5281/zenodo.21491997) 
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.1.0-blue.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📖 Overview
This repository contains the complete custom analytical pipeline, R scripts, and intermediate datasets for our study: *"A generalizable multi-pathway metabolic cell death signature predicts clinical outcomes and tumor microenvironmental states in bladder cancer"* (Currently under review/published in *Scientific Reports*).

Our study integrates transcriptomic data from multiple metabolic cell death (MCD) pathways to construct a robust, highly generalizable 8-gene prognostic signature for bladder cancer (BLCA). To ensure strict mathematical reproducibility and scientific restraint, all prognostic and immunotherapeutic inferences (including 1,000-iteration bootstrap resampling, cross-platform Z-score standardization, and *in silico* Immunophenoscore estimations) are fully open-sourced here.

## 📂 Repository Structure

The repository is organized into the following directories:

*   `01_Data/`: Contains formatted expression matrices and clinical datasets for the TCGA-BLCA training cohort and the GEO validation cohort (e.g., GSE13507). *(Note: Raw sequencing data should be downloaded directly from TCGA and GEO databases).*
*   `02_Scripts/`: R scripts for the complete analytical workflow.
    *   `Step1_DEG_Analysis.R`: Differential gene expression analysis and visualization (Volcano plots, Heatmaps).
    *   `Step2_Model_Construction.R`: LASSO Cox regression analysis to establish the 8-gene MCD signature.
    *   `Step3_Validation_and_Bootstrap.R`: Cross-cohort validation, Z-score standardization, and 1,000-iteration bootstrap resampling for C-index calculation.
    *   `Step4_Nomogram_Calibration.R`: Construction of the clinical-molecular nomogram and calibration curves.
    *   `Step5_TME_and_Immunotherapy.R`: ssGSEA for immune infiltration, Immunophenoscore (IPS) calculation, and drug sensitivity estimations.
    *   `Step6_Functional_Enrichment.R`: GO and KEGG pathway enrichment analyses (using `clusterProfiler`).
*   `03_Results/`: Output directory for generated figures, calibration curves, and supplementary tables.

## 🛠️ Prerequisites & Dependencies

The analytical pipeline was executed in **R (version 4.1.0 or higher)**. The following major R packages are required to run the scripts successfully:

```r
# Core Analysis & Machine Learning
install.packages(c("survival", "survminer", "glmnet", "rms", "pROC", "timeROC"))

# Data Manipulation & Visualization
install.packages(c("tidyverse", "ggplot2", "pheatmap", "ggpubr", "VennDiagram"))

# Bioinformatics & Enrichment (BiocManager)
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "GSVA", "limma"))

(Note: For drug sensitivity predictions, the pRRophetic package may require specific dependencies as outlined in its official documentation).

🚀 Usage Guide
To reproduce the findings of this study:

Clone this repository to your local machine:
git clone [https://github.com/sohnypaos-hue/metabolic-cell-death-related-gene-risk-model.git](https://github.com/sohnypaos-hue/metabolic-cell-death-related-gene-risk-model.git)
Open the R project file or set your working directory to the cloned repository.


