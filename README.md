# Integrated Multi-omics Analysis of Skeletal Muscle Ageing

![R](https://img.shields.io/badge/Language-R-276DC3)
![Data](https://img.shields.io/badge/Data-GTEx_v8-green)
![Analysis](https://img.shields.io/badge/Analysis-Multi--omics-orange)
![License](https://img.shields.io/badge/License-MIT-blue)

## Project Overview

This project investigates the biological mechanisms of human skeletal muscle ageing by integrating two distinct omics layers: **DNA methylation (DNAm)** and **Gene Expression (RNA-seq)**.

Using a meta-analysis of DNAm studies and RNA-seq. data from the **GTEx (v8) project**, we aimed to determine the extent to which epigenetic changes (methylation) couple with transcriptomic changes (expression) during the ageing process.

### Key Questions
1.  **Global Coupling:** Do methylation changes universally predict expression changes across the genome?
2.  **Directionality:** Is the relationship strictly canonical (Hypermethylation = Downregulation) or context-dependent?
3.  **Multi-omics Signatures:** Which genes show significant ageing effects in both molecular layers?

---

## Data Sources

| Modality | Source | Description |
| :--- | :--- | :--- |
| **DNA Methylation** | Meta-analysis | Aggregated EWAS summary statistics (BACONNED) for skeletal muscle. |
| **RNA-seq** | [GTEx Analysis v8](https://gtexportal.org/home/datasets) | "Muscle - Skeletal" tissue counts and phenotypes. |
| **Annotations** | Ensembl / HGNC | Mapped via `biomaRt` and custom annotation files. |

---

## Methodology

The analysis pipeline is implemented in R and follows these three stages:

### 1. DNA Methylation Processing
* **Input:** Probe-level summary statistics.
* **Aggregation:** Probes were mapped to genes; gene-level statistics were calculated (Mean Effect Size, Min P-value).
* **Filtering:** Focused on genes with valid symbols and sufficient coverage.

### 2. RNA-seq Analysis (`DESeq2`)
* **Preprocessing:** Downloaded GTEx gene count matrices and sample attributes.
* **Filtering:** Retained only "Muscle - Skeletal" samples and filtered low-count genes.
* **Modelling:** Modelled Age as an ordinal factor (`20-29` to `70-79`) to identify differentially expressed genes (DEGs).
* **Normalization:** VST (Variance Stabilizing Transformation) used for PCA and visualization.

### 3. Multi-omics Integration
* **Mapping:** Linked DNAm and RNA datasets using HGNC gene symbols.
* **Correlation:** Calculated Spearman correlations between DNAm Mean Effect Sizes and RNA Log2FoldChanges.
* **Categorization:** Grouped genes based on significance (FDR < 0.05) in one, both, or neither modality.

---

## Key Results

### 1. Global DNAm-RNA Coupling is Weak
Across all overlapping genes, the correlation between DNAm and RNA age effects was positive but weak (**ρ ≈ 0.124**). This suggests that global methylation changes explain only ~1.5% of the variance in age-related expression changes.

### 2. Categorisation of Ageing Genes
Among 22,791 overlapping genes:
* **41.4% (9,446)** were significant in **both** layers (Multi-omics hits).
* **28.3%** were significant in DNAm only.
* **10.0%** were significant in RNA only.

### 3. Directionality Patterns
Among the "Multi-omics hits" (significant in both), the direction of change varied:
* **40.0% Canonical (Inverse):** Hypermethylation/Downregulation or Hypomethylation/Upregulation.
* **60.0% Non-Canonical (Same-direction):** Hypermethylation/Upregulation or Hypomethylation/Downregulation.

### 4. Correlation by Category
The biological coupling is strongest in genes that are statistically significant drivers of ageing in both layers:
* **Both Significant:** ρ = 0.20
* **DNAm Only:** ρ = 0.08
* **RNA Only:** ρ = 0.03

> **Conclusion:** While global coupling is low, a specific subset of "multi-omics" genes drives the coordinated ageing response in skeletal muscle.

---

## Dependencies

To run this analysis, ensure you have **R (>= 4.0)** and the following packages installed:

```r
# CRAN
install.packages(c("tidyverse", "readxl", "metafor", "pheatmap", "ggplot2", "stringr", "readr"))

# Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "biomaRt"))