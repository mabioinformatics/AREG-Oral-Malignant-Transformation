# AREG-Oral-Malignant-Transformation
Overview
This repository contains reproducible R scripts investigating the association between amphiregulin (AREG) expression, glycolytic reprogramming, and malignant progression in oral potentially malignant disorders (OPMDs).

The analyses use publicly available transcriptomic and clinical data to examine whether elevated AREG expression is associated with increased glycolytic activity and subsequent oral cancer development.
Dataset
GSE26549 — Gene Expression Omnibus (GEO)
- Baseline OPMD biopsies: n = 62
- Data: Gene expression profiles and clinical follow-up
- Outcome: Oral cancer development
Dataset: GEO GSE26549
Analyses
Figure 4A. Association between AREG expression and glycolytic activity
Figure_4A_GSE26549.R

This analysis compares glycolysis scores between AREG-high and AREG-low baseline OPMDs.
- AREG expression: Mean expression of two transcript clusters (8095736 and 8095744).
- Glycolysis score: Mean of gene-wise standardized expression values (Z-scores) across a 10-gene signature.
- Group classification: Median AREG expression.
- Statistical analyses: Wilcoxon rank-sum test and Spearman correlation.
Figure 4B. AREG–glycolysis coupling and oral cancer-free survival
Figure_4B_GSE26549.R

This analysis evaluates the association between combined AREG-high/glycolysis-high status and oral cancer-free survival.
- AREG-high/Glycolysis-high: n = 17
- Other OPMDs: n = 45
- Follow-up: Administratively censored at 10 years
- Statistical analyses: Kaplan–Meier survival estimation, log-rank test, and Cox proportional hazards regression.
Glycolysis signature
The following 10 genes were used to calculate the glycolysis score:

SLC2A1, HK2, PFKP, ALDOA, ENO1, PKM, LDHA, PDK1, SLC16A1, and SLC16A3.

For each gene, expression values were standardized across baseline samples. The glycolysis score was calculated as the mean of the 10 standardized expression values for each sample.
Software Requirements
Analyses were performed using R with the following packages:
- GEOquery
- Biobase
- ggplot2
- survival
- survminer
Reproducibility
1. Download the R scripts from this repository.
2. Install the required R packages.
3. Run each script in RStudio.
4. The scripts automatically retrieve GSE26549 from GEO.
5. Statistical results are saved in the results/ directory.
6. Figures can be exported from the RStudio Plots panel.
The scripts use publicly available data and do not require access to restricted clinical information.
Data Availability
All gene expression and clinical data analyzed in this repository are publicly available through GEO under accession GSE26549.
