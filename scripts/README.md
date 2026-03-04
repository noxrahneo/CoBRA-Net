# Scripts Overview

This folder contains the first core scripts of the CoBRA-Net pipeline.

## 00_setup_environment.py

Purpose:
- Runs a sanity check before analysis starts.
- Validates configuration and required paths.
- Checks package availability.
- Writes a setup report to the logs folder.

Why it matters:
- Catches environment or path problems early.
- Makes runs more reproducible.

## 01_create_bigboss.py

Purpose:
- Builds the master sample inventory (`bigboss.csv`).
- Merges metadata from `SampleStats.txt` and Table EV4.
- Validates sample-name alignment between the two metadata sources.
- Adds matrix/barcodes filenames from raw data files.
- Outputs a clean, ordered table for downstream processing.

Why it matters:
- Creates one canonical metadata table used by all later analysis steps.

## 02_chopper.py

Purpose:
- Creates a filtered cohort file (`bigboss_chopped.csv`) from `bigboss.csv`.
- Applies thesis-scope filters (source, gender, condition).
- Removes samples outside the selected analysis scope.
- Prints a summary of kept rows.

Why it matters:
- Defines the exact cohort used in downstream analyses.

## Typical Run Order

From the repository root:

```bash
python3 scripts/00_setup_environment.py
python3 scripts/01_create_bigboss.py
python3 scripts/02_chopper.py
```

## Outputs

- `results/bigboss.csv` from `01_create_bigboss.py`
- `results/bigboss_chopped.csv` from `02_chopper.py`

## RtoPython Workflow Pointer

The current R-to-Python thesis workflow is in `scripts/RtoPython/`:

1. `01_pre_qc_from_r.py`
2. `02_filtering_cells.py`
3. `03_post_qc_from_r.py`
4. `04_per_sample_preprocess.py`
5. `05_validate_preprocess.py`
6. `06_integrate_per_condition.py`
7. `07_plot_integrated_results.py`
8. `08_plot_integrated_panels.py`
9. `09_annotate_clusters.py`
10. `10_annotate_clusters_minimal.py` (optional lightweight annotation)

See `scripts/RtoPython/README.md` for commands and outputs.

## Analysis Utilities

Non-translation analysis utilities live in `scripts/analysis/`:

- `04_marker_heatmaps.py`
- `05_composition_tests.py`
- `06_kegg_enrichment.py`
- `07_kegg_plots.py`
- `07_pseudobulk.py`
- `07c_pseudobulk_sanity_check.py`
- `07d_prepare_pre_correlation_pack.py`
- `07e_prepare_condition_logcpm_h5ad.py`
- `08a_compute_correlations.py`
- `08c_network_power_tom_prep.py`
- `08d_networkx_visualization.py`
- `slurm_run_08_network_array.sh`
- `08_pseudobulk_decoupler_downstream.py`
- `slurm_run_07_pseudobulk_array.sh`
- `slurm_run_08_network_build_and_viz.sh`
- `90_build_thesis_pack.py`
- `00_migrate_results_to_stages.py`
- `warehouse.py`

These utilities now read/write staged outputs under `results/stages/` by
default.

Pre-correlation preparation (recommended before network building):

```bash
python3 scripts/analysis/07_pseudobulk.py --condition all --threshold-mode auto
python3 scripts/analysis/07c_pseudobulk_sanity_check.py
python3 scripts/analysis/07d_prepare_pre_correlation_pack.py
python3 scripts/analysis/07e_prepare_condition_logcpm_h5ad.py
python3 scripts/analysis/08a_compute_correlations.py --condition all
python3 scripts/analysis/08c_network_power_tom_prep.py --condition all
python3 scripts/analysis/08d_networkx_visualization.py --condition all
```

Readable network visualization (interactive + cleaner static):

```bash
python3 scripts/analysis/08d_networkx_visualization.py \
	--condition all \
	--layout kamada \
	--max-edges 1500 \
	--min-weight 0.12 \
	--interactive-html \
	--interactive-max-edges 900 \
	--interactive-min-weight 0.18
```

Slurm array run for network prep + visualization:

```bash
sbatch scripts/analysis/slurm_run_08_network_array.sh
```

Optional Slurm overrides:
- `NETWORK_TYPE=signed|unsigned`
- `MAX_EDGES=2500`
- `MIN_WEIGHT=0.10`
- `POWERS=1,2,3,4,5,6,7,8,9,10,12,14,16,18,20`
- `INTERACTIVE_HTML=1`

Main network-ready inputs:
- `results/stages/07_network/pre_correlation/combined_logcpm.csv`
- `results/stages/07_network/pre_correlation/metadata_with_qc.csv`
- `results/stages/07_network/pre_correlation/include_mask_all.csv`
- `results/stages/07_network/pre_correlation/include_mask_exclude_flagged.csv`
- `results/stages/07_network/pre_correlation/per_condition/*_pseudobulk_logcpm.h5ad`

Migration:
- `python3 scripts/analysis/00_migrate_results_to_stages.py --mode copy`
- Optional dry run: `python3 scripts/analysis/00_migrate_results_to_stages.py --mode copy --dry-run`
