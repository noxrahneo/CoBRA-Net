# CoBRA-Net (Thesis Workspace)

This repository is the organized workspace for my thesis final version. It will evolve as the project progresses.

**Last updated**: February 16, 2026

## Current Structure

- `scripts/`: Main analysis workflow scripts (.py). These cover environment setup, data loading, preprocessing, QC analysis, and network building.
- `utils/`: Helper utilities used by the main scripts.

## Workflow Outline

1. Set up environment and load scRNA-seq data.
2. Preprocess data (basic scRNA-seq pipeline).
3. Initial data quality assessment (pre-QC) and visualization.
4. Perform QC and filtering (retain required data only).
5. Repeat quality assessment (post-QC) to show improvements and generate plots for the thesis.
6. Normalize data and select highly variable genes.
7. Annotate cell types.
8. Merge by cancer type and extract epithelial-only subsets:
	 - 1 healthy, 1 preneoplastic, 4 cancerous.
9. Compute correlations.
10. Build networks (individual and combined) to study cancer development.
11. Identify hub genes.
12. Downstream analyses: enrichment analysis, expression directionality,
		module comparison and therapeutic targeting to identify candidate drugs.

## Notes

- QC/preprocessing follows standard scRNA-seq best practices:
	https://www.sc-best-practices.org/introduction/analysis_tools.html
- Future additions may include garbage-collection helpers and expanded tooling.

## Current RtoPython Pipeline

Current active R-to-Python workflow (condition-aware):

1. `01_pre_qc_from_r.py`
2. `02_filtering_cells.py`
3. `03_post_qc_from_r.py`
4. `04_per_sample_preprocess.py`
5. `05_validate_preprocess.py`
6. `06_integrate_per_condition.py`
7. `07_plot_integrated_results.py`
8. `08_plot_integrated_panels.py`
9. `09_annotate_clusters.py`

From repo root:

```bash
python3 scripts/RtoPython/01_pre_qc_from_r.py --condition "ER+ tumor"
python3 scripts/RtoPython/02_filtering_cells.py --condition "ER+ tumor"
python3 scripts/RtoPython/03_post_qc_from_r.py --condition "ER+ tumor"
```

See `scripts/RtoPython/README.md` for details and staged outputs.

## Post-annotation analyses (all conditions)

From repo root:

```bash
python3 -m pip install --user gseapy statsmodels
python3 scripts/analysis/04_marker_heatmaps.py --condition all --group-col cell_type_annot
python3 scripts/analysis/05_composition_tests.py --condition all --group-col cell_type_annot --sample-col SampleName
python3 scripts/analysis/06_kegg_enrichment.py --condition all
```

Outputs:
- Marker heatmaps: `results/stages/04_annotation/<Condition>/figures/*marker_heatmap*.png`
- Composition test: `results/stages/05_composition/`
- KEGG enrichment: `results/stages/06_kegg/`
	- includes filtered interpretation tables: `*_kegg_interpretation.csv`

Provenance logs:
- Each stage writes a `warehouse.csv` in its output folder
	(`results/stages/04_annotation/`, `results/stages/05_composition/`,
	`results/stages/06_kegg/`, `results/stages/90_reports/thesis_pack/`).

Staged results:
- Non-destructive migration script: `python3 scripts/analysis/00_migrate_results_to_stages.py --mode copy`
- Staged layout index: `results/stages/README.md`

## Automation

This repo includes a `hooks/` folder to share git hook scripts. To enable them locally:

```bash
git config core.hooksPath hooks
```

Current hooks:

- `pre-commit`: updates the **Last updated** date in README.md.

I might add more hooks if needed