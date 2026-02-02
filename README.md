# TCGA-LIHC Multi-Omics Integration Pipeline

This repository contains a comprehensive R-based pipeline for the multi-omics analysis of Liver Hepatocellular Carcinoma (LIHC) data from TCGA. It covers the entire workflow from raw data acquisition to integrative clustering using the **MOVICS** framework.

##  Project Overview
The goal of this pipeline is to identify molecular subtypes of HCC by integrating five distinct layers of genomic and transcriptomic data. By finding the intersection of patients across all platforms, we ensure a high-resolution view of the tumor landscape.

### Integrated Data Types:
* **mRNA:** Protein-coding gene expression (Raw counts -> VST normalized).
* **lncRNA:** Long non-coding RNA expression.
* **DNA Methylation:** Epigenetic profiles (Illumina Human Methylation 450k).
* **Somatic Mutations:** Masked Somatic Mutation (MAF) converted to binary matrices.
* **CNA:** Copy Number Alterations (ASCAT2) sourced from GDC-Xena.

---

##  Workflow Architecture



### 1. Data Pre-processing
* **Normalization:** 
* **Cleaning:** 
* **ID Alignment:** Automated filtering to ensure all matrices contain the exact same 12-character TCGA patient barcodes.

### 2. Feature Selection (`getElites`)
To reduce noise, the pipeline selects "Elite" features:


### 3. Clustering & Subtyping
* **Consensus Clustering:** 
* **Evaluation:** Determines the optimal number of clusters ($k$) using:
    * Silhouette Scores
    * Gap Statistics
    * Consensus Matrices
### links:      
* https://github.com/xlucpu/MOVICS
* https://xlucpu.github.io/MOVICS/MOVICS-VIGNETTE.html
