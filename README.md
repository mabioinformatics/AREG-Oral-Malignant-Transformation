
AREG-Oral-Maligant-Transformation

# Amphiregulin and Glycolytic Reprogramming in Oral Malignant Transformation

## Overview

This repository provides reproducible R scripts for investigating the association between amphiregulin (AREG) expression, glycolytic reprogramming, and malignant progression in oral potentially malignant disorders (OPMDs).

The analyses use publicly available transcriptomic and clinical data to evaluate whether elevated AREG expression is associated with increased glycolytic activity and subsequent oral cancer development.

## Dataset

**GSE26549** — Gene Expression Omnibus (GEO)

Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26549

The analyses include 62 baseline OPMD biopsies with gene expression profiles and clinical follow-up information.

## Analyses

### Figure 4A. Association between AREG expression and glycolytic activity

**Script:** `Figure_4A_GSE26549.R`

This analysis compares glycolysis scores between AREG-high and AREG-low baseline OPMDs.

- AREG expression is calculated as the mean expression of two transcript clusters (8095736 and 8095744).
- AREG-high and AREG-low groups are defined using the median AREG expression.
- Glycolysis scores are calculated using a 10-gene signature.
- Statistical analyses include the Wilcoxon rank-sum test and Spearman correlation.

### Figure 4B. AREG–glycolysis coupling and oral cancer-free survival

**Script:** `Figure_4B_GSE26549.R`

This analysis evaluates the association between combined AREG-high/glycolysis-high status and oral cancer-free survival.

- AREG-high/Glycolysis-high: n = 17
- Other OPMDs: n = 45
- Follow-up is administratively censored at 10 years.
- Statistical analyses include Kaplan–Meier survival estimation, the log-rank test, and Cox proportional hazards regression.

## Glycolysis Signature

The glycolysis score is calculated using the following 10 genes:

SLC2A1, HK2, PFKP, ALDOA, ENO1, PKM, LDHA, PDK1, SLC16A1, and SLC16A3.

Expression values are standardized (Z-scores) for each gene across the 62 baseline samples. The glycolysis score for each sample is the mean of the 10 standardized expression values.

## Software Requirements

The analyses were developed in R using the following packages:

- GEOquery
- Biobase
- ggplot2
- survival
- survminer

## Reproducibility

1. Download the R scripts from this repository.
2. Install the required R packages.
3. Run each script in RStudio.
4. The scripts retrieve the GSE26549 dataset from GEO.
5. Statistical results and session information are saved in the `results/` directory.
6. Figures can be exported manually from the RStudio Plots panel.

An internet connection is required for the initial GEO data retrieval.

## Data Availability

All gene expression and clinical data used in these analyses are publicly available through GEO under accession GSE26549.

## Research Scope

These analyses explore associations between AREG expression, glycolytic activity, and malignant progression in OPMDs. The findings are observational and do not establish a causal relationship.

