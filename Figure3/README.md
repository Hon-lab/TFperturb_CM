# TF Perturbation Cardiomyocyte Time Course

## Overview

This repository contains analysis code for a TF perturbation cardiomyocyte time-course project, including cNMF gene program discovery, spatial MERFISH comparison, and embryo-stage comparison analyses.

The local workspace contains large AnnData objects, cNMF outputs, figures, and scheduler logs. The GitHub-ready repository is intended to track code and documentation only; generated data and large outputs remain on the analysis filesystem and are ignored by Git.

## Repository Structure

```text
.
├── cNMF_timecourse/       # Time-course cNMF runs, notebooks, scripts, configs, local outputs
├── spatial_comparison/    # MERFISH preprocessing, cNMF run, and spatial comparison notebooks
├── embryo_comparison/     # Embryo comparison notebooks and local outputs
├── Archive/               # Local metadata archived during repository cleanup
├── README.md
└── .gitignore
```

Within each analysis module, existing folder names are preserved to avoid breaking notebook-relative paths. Common folders include `code/` or `codes/` for notebooks/scripts, `data/` for local input AnnData/tables, `processed_data/` for generated tables, and `Plot/` or `figure/` for generated figures.

## Data

Large data files are not intended for Git tracking. Important local inputs include:

- `cNMF_timecourse/data/raw_count_D4.h5ad`
- `cNMF_timecourse/data/raw_count_D8.h5ad`
- `cNMF_timecourse/data/raw_count_D12.h5ad`
- `cNMF_timecourse/data/raw_count_D29.h5ad`
- `spatial_comparison/data/overall_merfish.h5ad`
- `spatial_comparison/data/overall_merfish_raw.h5ad`

Several notebooks also reference external HPC paths under `/project/shared/`, `/project/GCRB/`, and collaborator analysis directories. These paths are environment-specific dependencies and should be documented or parameterized before running the analysis on another system.

## Analysis Workflow

1. `cNMF_timecourse/preprocessing_data.ipynb` prepares raw count AnnData files for time-course cNMF.
2. `cNMF_timecourse/run_cNMF.sh` runs cNMF for D4, D8, D12, and D29 using the corresponding JSON config files.
3. `cNMF_timecourse/code/cNMF_analysis_*.ipynb` notebooks summarize cNMF programs and generate perturbation-level outputs.
4. `cNMF_timecourse/code/basic_stat_plot.ipynb`, `full_gp_stat.ipynb`, `gp_centrality_analysis.ipynb`, and `gp_interaction_analysis.ipynb` generate downstream summaries and figures.
5. `spatial_comparison/codes/Preprocessing_spatial.ipynb` prepares MERFISH data, and `spatial_comparison/run_cNMF.sh` runs spatial cNMF.
6. `spatial_comparison/codes/Gene_program_vis.ipynb` and `comapre_time_course_spatial.ipynb` compare time-course and MERFISH gene programs.
7. `embryo_comparison/code/` notebooks compare gene programs with embryo-stage datasets.

## Outputs

Generated outputs include cNMF result folders, processed CSV tables, AnnData intermediates, PDF figures, PNG figures, and Slurm logs. These are intentionally ignored by Git because the workspace is large and output files can be regenerated or stored separately.

## Notes

- The root `.gitignore` is configured for a lightweight GitHub repository with code and documentation only.
- Notebook path cleanup is not complete. Many notebooks still contain absolute HPC paths and should be refactored in a future reproducibility pass.
- The previous nested Git metadata from `embryo_comparison/.git` was archived under `Archive/embryo_comparison_git_metadata/` so `embryo_comparison/` can be treated as a normal folder in a root repository.
- No files were deleted during cleanup.
