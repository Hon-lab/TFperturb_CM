# TF Perturbation Gene Regulatory Network Analysis

## Overview

This repository contains Jupyter notebooks and generated outputs for analyzing transcription factor (TF) perturbation gene regulatory networks (GRNs) in a cardiac biology context. The notebooks combine pySpade perturbation hit tables, congenital heart disease (CHD) gene annotations, TF motif overlap tables, e-distance statistics, cNMF gene-program results, and graph analyses to build and summarize directed regulatory networks.

The main analysis products are:

- A directed GRN saved as `GRN_base.pickle`.
- Network topology, hierarchy, motif, PageRank, CHD, and gene-program summaries.
- PDF/PNG figures in `Plot/`.
- Small CSV summaries in `Data/` and `output_data/`.

Many primary inputs are not stored in this repository and are referenced by hard-coded absolute paths inside notebooks. Those private paths are intentionally not repeated here. To reproduce the full workflow, users will need to provide equivalent files and update the notebook paths.

## Repository Structure

| Path | Contents |
| --- | --- |
| `Archive/` | Earlier CHD-focused and cardiac-TF-focused exploratory notebooks. |
| `Code/` | Main notebooks for directed GRN construction, topology/statistics, gene-program analyses, and time-course analyses. |
| `Data/` | Small local CSV summaries and plotting data. |
| `output_data/` | Serialized layout objects and tabular outputs generated or copied by notebooks. |
| `Plot/` | Generated figures from the notebook workflow. |
| `GRN_base.pickle` | Local serialized NetworkX directed graph produced by `Code/Network_directed_analysis.ipynb`. |
| `plot.pdf`, `.png` | Standalone generated plot artifacts; their provenance is not clear from notebook references inspected here. |

## Complete Notebook Summary

| Notebook | Apparent role in workflow | Main inputs identified from code | Main outputs identified from active code | Notes and uncertainties |
| --- | --- | --- | --- | --- |
| `Code/Network_directed_analysis.ipynb` | Main directed GRN construction and CHD network analysis. Builds `GRN_base.pickle`, analyzes TF-TF hierarchy, integrates gnomAD LOEUF, computes reverse PageRank-like CHD regulator importance, and plots CHD regulator/target summaries. | External CHD gene table; external e-distance p-value table; external pySpade `filtered_df.csv`; external all-DE TF motif overlap table; gnomAD v4.1 constraint metrics TSV from Google Cloud. | `GRN_base.pickle`; `output_data/hierarchy_height.csv`; `Plot/example_network.pdf`; `Plot/hierarchy_DEG_LOEUF.pdf`; `Plot/hierarchy_simple.pdf`; `Plot/CHD_pagerank_importance.pdf`; `Plot/CHD_regulated_num.pdf`; `Plot/CHD_DE_num.pdf`. | This appears to be the core upstream notebook for the local graph and hierarchy output. It depends on external files and internet access for gnomAD unless cached. |
| `Code/Network_directed_stat.ipynb` | Network topology and motif statistics for `GRN_base.pickle`. Computes degree distributions, TF-TF motif census/randomization, and small-world style metrics. | Local `GRN_base.pickle`. | `Plot/Degree_stat_in.pdf`; `Plot/Degree_stat_out.pdf`; `Plot/motif_z_scores_barplot.pdf`; code writes `./data/motif_stats.csv`. | The repository contains `Data/motif_stats.csv` with capital `D`, while the notebook writes lowercase `data/motif_stats.csv`; this needs verification on case-sensitive filesystems. |
| `Code/Network_structure_analysis.ipynb` | Earlier/broader CHD GRN structure analysis. Builds regulatory graphs from CHD genes, e-distance hits, and pySpade hits, then analyzes clustering coefficients and path lengths for whole-network and TF-only views. | External CHD gene table; external e-distance p-value table; external pySpade `filtered_df.csv`. | `Plot/cluster_coef.pdf`; `Plot/ave_short_path.pdf`; `Plot/cluster_coef_TF.pdf`; `Plot/ave_short_path_TF.pdf`. | The notebook contains a commented CHD-TF GSEA save path. It appears exploratory and partly overlaps with `Network_directed_analysis.ipynb`. |
| `Code/Network_structure_gp.ipynb` | Gene-program GRN analysis for cNMF programs. Links perturbation effects to cNMF gene programs, builds per-program directed graphs, computes centrality/hub summaries, visualizes MYOCD-centered and top-target networks, and plots GO enrichment. | External pySpade `filtered_df.csv`; external cNMF regulation statistics pickle; external GO enrichment pickle; external top-300-genes-per-program CSV. | `output_data/gp_graph_pos.pickle`; `Plot/GRN_plot_with_legend.pdf`; `Plot/GRN_plot_total.pdf`; `Plot/GP_regulator_network_stat.pdf`; `Plot/GRN_MYOCD_centric.pdf`; `Plot/MYOCD_GP1_scatterplot.pdf`; `Plot/MYOCD_GP1_heatmap_Upstream_labels.pdf`; `Plot/GRN_plot_target.pdf`; `Plot/GP_GeneOntology_enrich.pdf`; `Plot/GP1_downstream_upstream_side_by_side.pdf`; `Plot/GP3_downstream_upstream_side_by_side.pdf`. | A notebook output shows a failed copy to `./output_data/top_300_genes_K_250_per_program_TF_full.csv` when run from `Code/`; the repository does contain that CSV under `output_data/` from the repository root. Path handling should be standardized. |
| `Code/Network_structure_time_course.ipynb` | Time-course extension of gene-program GRN analysis across D4, D8, D12, and D29 pySpade hit tables. Reuses gene-program results and saved graph positions to compare regulatory graphs over time, with MYOCD-focused summaries. | External static pySpade `filtered_df.csv`; external D4/D8/D12/D29 pySpade `filtered_df.csv` files; external cNMF regulation statistics pickle; external GO enrichment pickle; external top-300-genes CSV; local `output_data/gp_graph_pos.pickle`. | `Plot/GRN_plot_total_time_course.pdf`; `Plot/GRN_MYOCD_GP1_timecourse.pdf`; `Plot/GP_regulator_network_stat_time_course.pdf`; `Plot/GP_regulator_network_stat_time_course_top10.pdf`; `Plot/MYOCD_regulator_network_stat_time_course.pdf`; also writes/overwrites `Plot/GRN_plot_target.pdf` and `Plot/GP_GeneOntology_enrich.pdf`. | Some output filenames overlap with `Network_structure_gp.ipynb`, so execution order affects final files. |
| `Archive/20-5Disease_CHD-142genes.ipynb` | Archived CHD exploratory analysis. Performs CHD gene enrichment, intersects CHD genes with pySpade perturbation hits, builds CHD-TF and non-TF CHD subnetworks, and checks motif support. | External CHD gene table; external pySpade `filtered_df.csv`; external TSS BED; external non-TF CHD motif table; external `CHD_nonTF.gml` and position pickle. | `Plot/20-5GSEA_barchart_CHD_genes_v1.pdf` from active code. Several network figure save paths are present but commented. | Retained as exploratory/archive. Uses external graph/layout files not included here. |
| `Archive/24-6CardiacTFnetwork.ipynb` | Archived selected cardiac TF network analysis. Filters motif overlaps, builds selected cardiac-TF GRN, writes an external GML/layout, creates an interactive pyvis network, and analyzes regulated-gene overlap/enrichment. | External all-DE TF motif overlap table; external CHD gene table; external pySpade `filtered_df.csv`; external selected cardiac-TF GML/layout. | External GML and position pickle writes; `Plot/24-6GSEA_common_regulated_genes.pdf`; `Plot/24-6Venn5.pdf`; `Plot/24-6Venn6.pdf`; `Plot/24-6pseudovenn6.pdf`; pyvis output named `CardiacTF_gene_regulatory_network.html` in notebook output. | Main graph outputs were written outside the repository in the recorded notebook. The HTML file is not present in this repository. |
| `Archive/24-6CardiacTFnetwork_022525CT.ipynb` | Archived follow-up selected cardiac TF/CHD overlap analysis using the local base GRN. Compares sets of CHD genes, DEG genes, motif-supported targets, and local GRN structure. | External all-DE TF motif overlap table; external CHD gene table; local `GRN_base.pickle`; external pySpade `filtered_df.csv`. | `Plot/4gene_overlap.pdf`. | Appears to be a later review/reanalysis of the archived cardiac-TF work, but still depends on external primary tables. |

## Data Files

### Local Tabular and Serialized Files

| File | Size | Shape / object summary | Apparent columns or contents | Apparent role |
| --- | ---: | --- | --- | --- |
| `Data/gene_estat_emb_tsne_data.csv` | 6.8 KB | 209 rows x 5 columns | `sgRNA`, `x`, `y`, `cluster`, plus an index column | Appears to contain two-dimensional embedding coordinates and cluster IDs for sgRNA/gene perturbations. I did not find a notebook reference to this file, so its exact workflow role needs verification. |
| `Data/motif_stats.csv` | 1.2 KB | 16 rows x 8 columns | motif triad code, real count, random mean/std, z-score, over/under empirical p-values | Motif census summary corresponding to `Network_directed_stat.ipynb`. The notebook writes lowercase `data/motif_stats.csv`, while this repository stores `Data/motif_stats.csv`. |
| `output_data/hierarchy_height.csv` | 3.1 KB | 83 rows x 7 columns | gene/TF index, `hierarchy_height`, in/out/total degree, CHD-TF flag, LOEUF | Final/intermediate output from `Network_directed_analysis.ipynb` for TF hierarchy and constraint analyses. |
| `output_data/top_300_genes_K_250_per_program_TF_full.csv` | 584.7 KB | 300 rows x 250 columns | columns `1` to `250`, each listing ranked genes for a cNMF gene program | Local copy of a cNMF top-gene table used by the gene-program notebooks. The notebooks currently read an external path, not this relative copy. |
| `GRN_base.pickle` | 1.3 MB | NetworkX `DiGraph`; 8,438 nodes and 25,621 edges | node attributes include `type`, `color`, `subset`, `motif_info`; edge attributes include `regulation`, `is_motif` | Core local GRN object generated by `Network_directed_analysis.ipynb` and reused by downstream notebooks. |
| `output_data/gp_graph_pos.pickle` | 178.5 KB | Python dictionary with 232 graph-layout entries | keys such as `Usage_1_total`; values are node-coordinate dictionaries | Saved graph layouts from `Network_structure_gp.ipynb`, reused by `Network_structure_time_course.ipynb`. |

### External Inputs Referenced by Notebooks

These inputs are important to the workflow but are not included in this repository. The table lists logical filenames and inspected structure only, not private absolute paths.

| Logical input | Inspected size / shape | Main columns or object contents | Used by | Apparent role |
| --- | --- | --- | --- | --- |
| `chdgene_table.csv` | 11.7 KB; 142 rows x 6 columns | `Gene`, CHD classification, extracardiac phenotype, inheritance mode, ranking, supporting references | CHD and directed-network notebooks | CHD disease gene annotation table. |
| Static pySpade `filtered_df.csv` | 5.3 MB; 27,578 rows x 18 columns | hit index, `gene_names`, genomic coordinates, perturbation `region`, cell/bin counts, hypergeometric p-value score, `fc`, `Significance_score`, empirical p-value, perturb/background CPM | Most notebooks | Primary perturbation hit table linking perturbed regions/TFs to regulated genes. |
| `p_val_gRNA_sig_full_expressing_v3.csv` | 108.9 KB; 211 rows x 50 columns | repeated distance/p-value columns, p-value and distance summary statistics, `gene_name`, sgRNA/cell counts | `Network_directed_analysis.ipynb`, `Network_structure_analysis.ipynb` | e-distance significance table used to identify valid TF perturbations. |
| `All_DEgenes_TFmotif.txt` | 1.79 GB; 23,141,698 rows x 12 no-header columns | columns assigned in notebooks as regulatory element coordinates, target `gene`, motif ID/position/strand, `TF`, overlap length | Directed and cardiac-TF notebooks | Motif-overlap support for TF-target regulatory edges across all differentially expressed genes. |
| `NonTFCHD_gene_TFmotif.txt` | 23.1 MB; 301,985 rows x 12 no-header columns | same motif-overlap schema as above | `Archive/20-5Disease_CHD-142genes.ipynb` | Motif-overlap support for non-TF CHD gene network edges. |
| Time-course pySpade `D4/D8/D12/D29 filtered_df.csv` | D4: 38,466 rows; D8: 69,172 rows; D12: 69,076 rows; D29: 79,771 rows; each 18 columns | same schema as static pySpade `filtered_df.csv` | `Network_structure_time_course.ipynb` | Time-point-specific perturbation hit tables for longitudinal gene-program GRNs. |
| `compare_score_pert_all_k_all_targets_clear_gRNA_12_3_24.pkl` | 174.3 MB; dictionary with 32 keys | keys such as `usage_norm_k_5`; values are DataFrames with `Gene`, score column, statistic, p-value, log2 fold change, BH-corrected p-value | Gene-program notebooks | Perturbation effects on cNMF usage/gene-program scores. |
| `TF_full_GO_all_enriched_results.pkl` | 463.1 MB; dictionary with 32 keys | values are enrichment DataFrames with gene set, term, overlap, p-values, odds ratio, combined score, genes, source | Gene-program notebooks | GO/enrichment annotations for cNMF gene programs. |

## How the Data Fit Together

The workflow appears to start from TF perturbation hit calls in pySpade `filtered_df.csv`, where each row links a perturbation region to a regulated gene and records effect size/significance. Annotation dictionaries inside notebooks map perturbation regions back to genes/TFs.

CHD analyses add `chdgene_table.csv` to identify disease genes and disease-associated TFs among perturbation targets. The directed-network notebooks combine pySpade hits, e-distance significant perturbations, and motif-overlap support from `All_DEgenes_TFmotif.txt` to construct `GRN_base.pickle`.

`GRN_base.pickle` is then used for graph topology, motif census, hierarchy, LOEUF, and CHD regulator analyses. `output_data/hierarchy_height.csv` is a derived table from this graph, and `Data/motif_stats.csv` is a derived motif-randomization summary.

Gene-program notebooks connect perturbation effects to cNMF programs using the external regulation-statistics pickle and top-gene table. They build program-specific TF-to-gene-program or TF-to-target graphs, save graph layouts to `output_data/gp_graph_pos.pickle`, and generate summary PDFs. The time-course notebook applies the same logic across D4, D8, D12, and D29 pySpade tables.

## Workflow

An inferred execution order is:

1. Prepare external inputs: CHD gene annotations, pySpade hit tables, e-distance statistics, TF motif overlaps, cNMF program statistics, cNMF top genes, and GO enrichment results.
2. Run `Code/Network_directed_analysis.ipynb` to construct `GRN_base.pickle`, write `output_data/hierarchy_height.csv`, and generate CHD/directed-network plots.
3. Run `Code/Network_directed_stat.ipynb` to compute graph topology and motif statistics from `GRN_base.pickle`.
4. Optionally run `Code/Network_structure_analysis.ipynb` for exploratory whole-network and TF-only clustering/path-length analyses.
5. Run `Code/Network_structure_gp.ipynb` for gene-program GRNs and graph layouts.
6. Run `Code/Network_structure_time_course.ipynb` for time-course gene-program analyses using the saved graph layouts.
7. Use `Archive/` notebooks only when reproducing historical exploratory analyses.

This order is inferred from reads/writes and should be verified before rerunning, because several notebooks contain hard-coded paths, relative-path assumptions, commented save statements, and overlapping output filenames.

## Input Data

The repository does not include most raw or primary analysis inputs. Required external inputs include:

- CHD gene annotation table.
- Static and time-course pySpade `filtered_df.csv` perturbation hit tables.
- e-distance p-value/significance table.
- All-DE and CHD-specific TF motif overlap tables.
- cNMF perturbation score/statistics pickle.
- cNMF GO/enrichment pickle.
- Top 300 genes per cNMF program table.
- gnomAD constraint metrics TSV, downloaded by URL in `Code/Network_directed_analysis.ipynb`.

For reproducibility, copy these inputs into a documented local `data/` directory or configure a path file, then replace the hard-coded notebook paths.

## Output Files

Important generated outputs include:

- `GRN_base.pickle`: directed GRN used by downstream analyses.
- `output_data/hierarchy_height.csv`: TF hierarchy/degree/LOEUF table.
- `output_data/gp_graph_pos.pickle`: saved layouts for gene-program graphs.
- `Data/motif_stats.csv`: motif census/randomization statistics.
- `Plot/*.pdf` and `Plot/*.png`: network, enrichment, topology, CHD, gene-program, and time-course figures.

Some notebooks write the same output filenames, including `Plot/GRN_plot_target.pdf`, `Plot/GP_GeneOntology_enrich.pdf`, and `Plot/MYOCD_regulator_network_stat_time_course.pdf`. Rerunning notebooks can overwrite existing figures.

## Installation and Environment

The recorded notebook outputs show Python 3.12.2 for most current notebooks and Python 3.7 for at least one archived notebook. Package versions recorded in several current notebooks include:

- `numpy` 1.26.4
- `pandas` 2.2.1
- `scipy` 1.15.3
- `matplotlib` 3.8.3

Imports used across notebooks include:

- `networkx`
- `gseapy`
- `seaborn`
- `tqdm`
- `joblib`
- `pyvis`
- `venn`
- `matplotlib`
- `pandas`
- `numpy`
- `scipy`

No environment file is currently provided. A practical starting point is to create a conda or mamba environment with the packages above plus Jupyter.

## Usage

1. Create and activate a Python/Jupyter environment with the required packages.
2. Place required external inputs in accessible locations.
3. Update notebook path variables to point to those inputs. Avoid committing private absolute paths.
4. Run notebooks from the repository root unless a notebook explicitly assumes another working directory. Several relative paths differ between `Code/` and the repository root.
5. Start with `Code/Network_directed_analysis.ipynb` if rebuilding the core GRN.
6. Run downstream statistics and gene-program notebooks after confirming `GRN_base.pickle` and `output_data/gp_graph_pos.pickle` are present as needed.

## Reproducibility Notes

- The workflow is notebook-driven and not yet packaged as scripts or a pipeline.
- Hard-coded external paths must be replaced for portability.
- Large motif tables should be processed by streaming/chunked methods where possible; the all-DE motif table is approximately 1.79 GB.
- The archived notebooks include useful provenance but may depend on older environments and external graph/layout files.
- Some save paths are commented, and some files in `Plot/` may come from prior runs rather than active code in the current notebooks.
- The graph contains edge attribute values spelled `up_regualte` and `down_regualte` in the serialized object; downstream code should preserve or normalize these carefully.
- Several notebooks use randomization/layout algorithms. Set random seeds consistently before regenerating figures intended for exact reproduction.
- Notebook execution order, working directory, and output overwrites should be documented before publication or handoff.

## Citation

No publication citation or formal software citation file was found in the repository. If this repository supports a manuscript, cite the associated study and any external tools/resources used, including pySpade, NetworkX, gseapy/Enrichr, cNMF, gnomAD, and the CHD gene source.

## Contact

Chikara Takeuchi (chikara.takeuchi@utsouthwestern.edu) - UT Southwestern Medical Center
