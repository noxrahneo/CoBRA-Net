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

See `scripts/RtoPython/README.md` for commands and outputs.

## Analysis Utilities

Non-translation analysis utilities live in `scripts/analysis/`:

- `04_marker_heatmaps.py`
- `05_composition_tests.py`
- `06_kegg_enrichment.py`
- `90_build_thesis_pack.py`
- `00_migrate_results_to_stages.py`
- `warehouse.py`

These utilities now read/write staged outputs under `results/stages/` by
default.

Migration:
- `python3 scripts/analysis/00_migrate_results_to_stages.py --mode copy`
- Optional dry run: `python3 scripts/analysis/00_migrate_results_to_stages.py --mode copy --dry-run`
