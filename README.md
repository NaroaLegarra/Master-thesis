---
title: "README"
author: "Naroa Legarra Marcos"
date: "2022-08-23"
output: html_document
---

# Predicting the cellular composition and molecular makeup of the neurovascular unit using machine learning methods

This repository contains the code used in the course of the master's thesis. 

## Code

### Dimensionality reductions

- *dimensionality_reductions.py* : Script to make PCA+TSNE in sc, pseudobulk and pseudobulk replicas with the CUN dataset, Granja dataset and Processed Granja dataset. 

### DEA analysis

- *dea_analysis_data_processing.R* :Process the sc RNAseq and FACS-bulk RNAseq data from the CUN dataset data for further analysis

- *dea_analysis_functions.R* : Functions used to process and perform the analysis of sc RNAseq and FACS-bulk RNAseq data

- *dea_analysis_main.R* : Performs the differential gene expression analysis between cell types in FACS-bulk RNAseq and pseudobulk RNAseq in the CUN datasets and compares results

- *dea_analysis_plots.R* : Generates plots for marker gene, enrichment analysis and DEA comparison

### CIBERSORTx

- *cibersortx_processing.R* : Process the rds files to save the psuedobulk matrices in tsv format for cibersortx

- *cibersortx_top_genes_selection.R* : Selects the top differentially expressed genes for each type of cell vs the other two using FACS-bulk RNAseq data

- *cibersortx_heatmaps.py* : Script to make heatmaps comparing reference and CIBERSORTx predicted results
