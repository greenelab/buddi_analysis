# BuDDI Analysis

This repository contains all the code to reproduce the analyses in [**BuDDI: Bulk Deconvolution with Domain Invariance to predict cell-type-specific perturbations from bulk**](https://www.biorxiv.org/content/10.1101/2023.07.20.549951v1)

## Directory Structure
Data for each analysis can be downloaded from figshare under the DOI: 10.6084/m9.figshare.23721336.
The expected directory structure is:
- buddi_analysis
    - sc_preprocessing
    - kang_analysis
    - liver_analysis
    - synovium_analysis
    - bayesprism_scripts
    - data (from figshare)
    - results (from figshare)

## Contents
A tutorial on how to run a BuDDI analysis is provided in [tutorial folder](https://github.com/greenelab/buddi_analysis/tree/main/tutorial).
Data to run the tutorial are also provided in the same folder.

The analysis folders contain the following information:
- [sc_preprocessing](https://github.com/greenelab/buddi_analysis/tree/main/sc_preprocessing): single-cell QC scripts and pseudobulk generation code
- [kang_analysis](https://github.com/greenelab/buddi_analysis/tree/main/kang_analysis): code to generate figures 1-3, supp. fig 1,5,6, supp. table 1
    - This includes the code for CVAE and PCA comparisons in the folder [kang_comparators](https://github.com/greenelab/buddi_analysis/tree/main/kang_analysis/kang_comparators)
- [liver_analysis](https://github.com/greenelab/buddi_analysis/tree/main/liver_analysis): code to generate figure 4, supp. fig 2,3,6, supp. table 2,3
    - This includes the code for CVAE and PCA comparisons in the folder [sex_prediction_comparators](https://github.com/greenelab/buddi_analysis/tree/main/liver_analysis/sex_prediction_comparators)
- [synovium_analysis](https://github.com/greenelab/buddi_analysis/tree/main/synovium_analysis): code to generate figure 5, supp. fig 7, supp. table 4
- [bayesprism_scripts](https://github.com/greenelab/buddi_analysis/tree/main/bayesprism_scripts): R and shell scripts to run [BayesPrism](https://doi.org/10.1038/s43018-022-00356-3)


## Running Analyses
To use these notebooks, it is assumed that you have installed BuDDI in a virtualenv and have activated this venv before running the notebooks.
Installation instructions for BuDDI are here: https://github.com/greenelab/buddi/tree/main#buddi

