# R to Python Conversions

This folder contains Python conversions of selected scripts from `HumanBreast10X-main/RCode`.

## Why this folder exists

- Keep direct R→Python translations in one place.
- Preserve the logic of the original publication code.
- Make conversions easy to review, compare and cite in thesis methods.

## Current conversions

### `pre_qc_from_r.py`

Purpose:
- Runs pre-filter QC directly on raw matrices from `results/bigboss_chopped.csv`.
- Computes sample-level and per-cell QC metrics before filtering.

Outputs:
- `results/r_to_python/qc_pre/pre_qc_summary.csv`
- `results/r_to_python/qc_pre/pre_qc_per_cell.parquet` (or CSV fallback)
- `results/r_to_python/qc_pre/pre_qc_boxplots.png`
- `results/r_to_python/qc_pre/pre_qc_violins.png`

Run:
- `python scripts/RtoPython/pre_qc_from_r.py`
- `python scripts/RtoPython/pre_qc_from_r.py --condition "ER+ tumor"`

### `filtering_cells.py`

Source R script:
- `data/HumanBreast10X-main/RCode/NormTotal.R`

What is implemented:
- Minimal **cohort filtering** workflow using `results/bigboss_chopped.csv`.
- Loads each selected sample from raw matrix + barcodes + shared features file.
- Computes per-cell QC metrics required by R logic:
	- library size (`lib_size`)
	- detected genes (`n_genes`)
	- mitochondrial fraction (`percent_mito`)
- Applies sample-specific cell thresholds from metadata:
	- `percent_mito < Mito`
	- `GeneLower < n_genes < GeneUpper`
	- `lib_size < LibSize`
- Applies basic gene filtering:
	- keep genes detected in at least 1% of cells
	- remove empty/invalid symbols
	- remove duplicated gene symbols
- Saves new filtered matrices as `.h5ad` files (one per sample).
- Works on all conditions by default; optional `--condition` lets you run one condition only.
- Writes a summary CSV with before/after cell and gene counts.

Outputs (default):
- `results/filtered_samples/<Condition>/<SampleName>_filtered.h5ad`
- `results/filtered_samples/filtering_summary.csv`

Run:
- `python scripts/RtoPython/filtering_cells.py`
- `python scripts/RtoPython/filtering_cells.py --condition Normal`
- `python scripts/RtoPython/filtering_cells.py --list-conditions`

Important clarification:
- This script performs **actual filtering** and writes new filtered outputs.
- It is still a scoped conversion focused on filtering logic, not a full reproduction of every downstream Seurat step in `NormTotal.R`.

### `post_qc_from_r.py`

Purpose:
- Runs post-filter QC on filtered `.h5ad` files from `results/filtered_samples`.
- Uses the same metric definitions as pre-QC.
- Compares pre/post QC when pre summary is available.

Outputs:
- `results/r_to_python/qc_post/post_qc_summary.csv`
- `results/r_to_python/qc_post/post_qc_per_cell.parquet` (or CSV fallback)
- `results/r_to_python/qc_post/post_qc_boxplots.png`
- `results/r_to_python/qc_post/post_qc_violins.png`
- `results/r_to_python/qc_post/pre_post_qc_comparison.csv` (if pre summary exists)
- `results/r_to_python/qc_post/pre_post_qc_scatter.png` (if pre summary exists)
- `results/r_to_python/qc_post/pre_post_qc_violins.png` (if pre per-cell exists)

Run:
- `python scripts/RtoPython/post_qc_from_r.py`
- `python scripts/RtoPython/post_qc_from_r.py --condition "ER+ tumor"`

## Notes

- This conversion reproduces the QC visualization logic; it does not run full Seurat workflows.
- Additional R scripts (annotation/clustering/group-specific scripts) can be converted here incrementally.
