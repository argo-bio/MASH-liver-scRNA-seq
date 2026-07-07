# Single-Cell RNA Sequencing Analysis Toolkit

## Overview
This repository contains the R codebase for analyzing single-cell RNA sequencing (scRNA-seq) data derived from mouse models of Metabolic Dysfunction-Associated Steatohepatitis (MASH). The project investigates the role of the **Wnt/β-catenin** signaling pathway in hepatic macrophage immunometabloism reprograming and disease progression.

The workflow is divided into three main components: benchmarking using public datasets, primary analysis of in-house generated data (Ctnnb1 CKO vs. ACT), and a library of reusable R functions for standardized single-cell workflows.


## Directory Structure
├── public-MASH-liver-scRNA-seq/     # Analysis of public reference datasets
│   ├── step1_public_scRNA-seq-anno...   # Raw data processing & cell annotation
│   └── step2_enrichment_macrophag...    # Macrophage subset enrichment analysis
│
├── Ctnnb1-CKO-ACT-MASH-liver/       # Core analysis (In-house data)
│   ├── step1_total_hepatic_landscape.R # Hepatic cellular landscape construction
│   ├── step2_immune_cell_porotion_...   # Immune cell proportion statistics
│   ├── step3_enrichment_analysis.R     # DEG & Pathway enrichment analysis
│   └── step4_cellchat_mac-hep_LESCs.R  # Cell-cell communication (CellChat)
│
├── sc_workflow_QingLab/             # Custom R function library
│   ├── run_harmony.R                # Batch correction via Harmony
│   ├── run_cellchat.R               # Automated CellChat pipeline
│   ├── run_marker.R                 # Differential expression & marker identification
│   └── run_Vln.R                    # Visualization utilities (Violin plots, etc.)
│
└── README.md                        # This file

## Environment & Dependencies
All analyses were performed in R (>= 4.0).

## Data Availability

Due to size constraints, raw sequencing data (FASTQ) are not included in this repository. Processed count matrices (e.g., `.mtx` or `.rds`) required to run these scripts are available from the corresponding author upon reasonable request.

## Contact

For questions regarding the code, please open an issue on GitHub or contact the maintainer.
