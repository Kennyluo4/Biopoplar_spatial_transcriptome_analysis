# Biopoplar_spatial_transcriptome_analysis

Scripts and computational workflows for spatial transcriptomics analysis of poplar 717 spatial transcriptome.

## Overview

This repository contains scripts and analysis pipelines used to process, quantify, and analyze spatially resolved transcriptomic data. The workflows encompass raw data processing, spatial gene expression profiling, clustering, and downstream functional annotations.

## Requirements & Dependencies

* **R** ($\ge 4.0.0$)
  * `Seurat` / `SeuratObject`
  * `ggplot2`
  * `dplyr`
  * `patchwork`
* **Python** ($\ge 3.8$)
  * `scanpy`
  * `spatialdata`
  * `pandas`
  * `numpy`
  * `matplotlib` / `seaborn`

## Repository Structure

```text
├── R/
│   ├── common_helpers.R
│   └── submission_setup.R
├── resources/
├── 01_shoot_apex.R
├── 02_axillary_bud.R
├── 03_stem.R
├── 04_petiole_cross.R
├── 05_petiole_longitudinal.R
├── 06_generate_figures.R
└── README.md