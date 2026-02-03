# Differential Gene Expression Analysis (limma-voom)

## Overview
This repository contains a complete **limma-voom differential gene expression workflow** for bulk RNA-seq data derived from a **rat feeding study** comparing two dietary conditions.

The analysis focuses on identifying transcriptional differences between rats fed **corn oil** and those fed **flaxseed oil**, using established Bioconductor methods.

---

## Experimental Design
- **Organism:** Rat  
- **Data type:** Bulk RNA-seq (gene-level counts)
- **Groups:**
  - **LIC** – rats fed **corn oil**
  - **LIF** – rats fed **flaxseed oil**
- **Replicates:** 5 biological replicates per group

The comparison performed is:

> **Flaxseed-fed (LIF) vs Corn oil-fed (LIC)**

---

## Analysis Workflow
The analysis was conducted using the **limma-voom** framework, which combines the strengths of `edgeR` and `limma`:

1. Import raw gene-level count matrix
2. Align count data with sample metadata
3. Normalize counts using **TMM normalization**
4. Apply **voom transformation** to model the mean–variance relationship
5. Fit linear models and apply empirical Bayes moderation
6. Identify differentially expressed genes
7. Generate standard RNA-seq visualizations

---

## Outputs
The workflow produces the following outputs:

- **Differential expression results**
  - `results/limma_voom_results.csv`
- **Visualizations**
  - Volcano plot
  - MA plot
  - Heatmap of top 30 differentially expressed genes

These plots are saved in the `figures/` directory.

---

## Repository Structure
