# Gastric Myeloid Transcriptomics

## Overview

This repository contains R analysis code supporting the manuscript *Integrated single-cell and bulk transcriptomic analysis identifies macrophage reprogramming-related prognostic markers in gastric carcinogenesis*.

The study integrates single-cell and bulk transcriptomic analyses to investigate immune microenvironment remodeling, macrophage-associated cellular heterogeneity, transcriptional regulatory programs, and survival-associated candidate genes in gastric cancer.

## Analysis workflow

The main analysis file is `PLOS_ONE_Methods_Fig1_Fig6.R`.

The workflow includes:

- M00–M03: Software environment, data inspection, single-cell annotation, and cell-type composition.
- M04–M04C: Pseudobulk differential expression, Hallmark GSEA, and tissue-sample-level program scores.
- M05: CellChat communication network analysis.
- M06–M07: Macrophage reclustering, doublet detection, ambient RNA assessment, and quality-control sensitivity analyses.
- M08: CIBERSORT results and cross-platform comparisons.
- M09: Slingshot trajectory inference and root-selection sensitivity analyses.
- M10–M11: SCENIC AUCell activity, transcription factor networks, and functional enrichment.
- M12: Cox regression and external cohort survival analyses.
- M13: Reproducibility and input verification.

## Data sources

Public datasets include GSE134520, GSE234129, GSE191275, GSE62254, GSE84437, and TCGA-STAD data obtained through UCSC Xena.

Original large expression matrices and restricted third-party sequencing data are not redistributed in this repository.

## Reproducibility status

This repository contains a consolidated workflow incorporating existing analysis routines, saved-result verification, and reconstructed computational modules.

Selected modules have been tested against locally available datasets and intermediate analysis objects. Certain computational modules require additional verified inputs or comparison with the final manuscript results.

In particular, modules M04, M04C, M11, and M12 have unresolved input-verification requirements, while M06 and M09 require final result comparisons.

Successful execution of an individual verification module does not necessarily indicate full end-to-end reproduction of its associated figure.

The workflow requires configuration of local input paths and availability of the corresponding analysis data.

## Software

Analyses were conducted in R using packages including Seurat, Harmony, edgeR, limma, CellChat, scDblFinder, celda/DecontX, Slingshot, SCENIC, AUCell, clusterProfiler, and survival.

Detailed software versions should be obtained from the accompanying R session information.

## Code availability

The scripts are provided to support methodological transparency and reproducibility. Additional verification records and analysis inputs will be incorporated as they become available.
