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

### `validate_preprocess.py`

Purpose:
- Validate per-sample preprocessing outputs from `per_sample_preprocess.py`.
- Compare filtered vs preprocessed files sample-by-sample.
- Check preprocessing artifacts and consistency:
	- `counts` layer exists
	- `highly_variable` column exists
	- HVG count passes selected mode (`strict` or `relaxed`)
	- no non-finite values in expression matrix
	- cell/gene dimensions are preserved

Outputs:
- Per-condition CSV:
	- `results/r_to_python/preprocess_validation/preprocess_validation_<Condition>.csv`
- Run-level summary JSON:
	- `results/r_to_python/preprocess_validation/preprocess_validation_<Condition>.json`
- If `--condition all`, also writes:
	- `results/r_to_python/preprocess_validation/preprocess_validation_all.csv`
	- `results/r_to_python/preprocess_validation/preprocess_validation_all.json`

Run:
- List available condition folders:
	- `python scripts/RtoPython/validate_preprocess.py --list-conditions`
- Validate one condition:
	- `python scripts/RtoPython/validate_preprocess.py --condition Normal`
- Validate one condition (strict HVG mode):
	- `python scripts/RtoPython/validate_preprocess.py --condition Normal --hvg-mode strict --expected-hvg 3000`
- Validate one condition (relaxed HVG mode):
	- `python scripts/RtoPython/validate_preprocess.py --condition Normal --hvg-mode relaxed --hvg-tolerance 0.05 --expected-hvg 3000`
- Validate all conditions:
	- `python scripts/RtoPython/validate_preprocess.py --condition all --hvg-mode relaxed`

### `per_sample_preprocess.py`

Purpose:
- Run per-sample preprocessing before integration.
- Seurat-style step mapping: NormalizeData -> FindVariableFeatures -> ScaleData.
- Support one condition or all conditions from condition folders.

Outputs:
- Preprocessed sample files:
	- `results/preprocessed_samples/<Condition>/<SampleName>_preprocessed.h5ad`
- Summary CSV:
	- one condition: `results/preprocessed_samples/<Condition>/preprocess_summary.csv`
	- all conditions: `results/preprocessed_samples/preprocess_summary.csv`

Run:
- List available condition folders:
	- `python scripts/RtoPython/per_sample_preprocess.py --list-conditions`
- Run one condition:
	- `python scripts/RtoPython/per_sample_preprocess.py --condition Normal`
- Run one condition by label (safe folder mapping supported):
	- `python scripts/RtoPython/per_sample_preprocess.py --condition "ER+ tumor"`
- Run all conditions:
	- `python scripts/RtoPython/per_sample_preprocess.py --condition all`

### `integrate_per_condition.py`

Purpose:
- Integrate preprocessed samples per condition into a joint condition object.
- Seurat-style mapping:
	- FindIntegrationAnchors / IntegrateData (Python equivalent)
	- joint PCA / neighbors / UMAP / Leiden
- Supports one condition or all conditions.

Outputs:
- Integrated condition object:
	- `results/integrated_samples/<Condition>/<Condition>_integrated.h5ad`
- Per-cell metadata table:
	- `results/integrated_samples/<Condition>/<Condition>_cell_metadata.csv`
- UMAP coordinates table:
	- `results/integrated_samples/<Condition>/<Condition>_umap.csv`
- Integration summary:
	- one condition: `results/integrated_samples/integration_summary_<Condition>.csv`
	- all conditions: `results/integrated_samples/integration_summary.csv`

Run:
- List available condition folders:
	- `python scripts/RtoPython/integrate_per_condition.py --list-conditions`
- Integrate one condition:
	- `python scripts/RtoPython/integrate_per_condition.py --condition Normal`
- Integrate one condition by label (safe folder mapping supported):
	- `python scripts/RtoPython/integrate_per_condition.py --condition "ER+ tumor"`
- Integrate all conditions:
	- `python scripts/RtoPython/integrate_per_condition.py --condition all`

### `plot_integrated_results.py`

Purpose:
- Create post-integration visualization plots per condition.
- Save UMAP plots by cluster (`leiden`) and by sample (`SampleName`).
- Optional t-SNE plotting for article-style embedding comparison.

Inputs:
- `results/integrated_samples/<Condition>/<Condition>_integrated.h5ad`

Outputs:
- `results/integrated_samples/<Condition>/figures/<Condition>_umap_leiden.png`
- `results/integrated_samples/<Condition>/figures/<Condition>_umap_sample.png`
- Optional (`--compute-tsne`):
	- `results/integrated_samples/<Condition>/figures/<Condition>_tsne_leiden.png`
	- `results/integrated_samples/<Condition>/figures/<Condition>_tsne_sample.png`

Run:
- List available condition folders:
	- `python scripts/RtoPython/plot_integrated_results.py --list-conditions`
- Plot one condition:
	- `python scripts/RtoPython/plot_integrated_results.py --condition Normal`
- Plot all conditions:
	- `python scripts/RtoPython/plot_integrated_results.py --condition all`
- Plot all conditions with t-SNE:
	- `python scripts/RtoPython/plot_integrated_results.py --condition all --compute-tsne`
- If legends are crowded/cut, set legend options:
	- `python scripts/RtoPython/plot_integrated_results.py --condition Normal --legend-loc "right margin" --legend-fontsize 7`

### `plot_integrated_panels.py`

Purpose:
- Build publication-style multi-condition UMAP panel figures.
- Use fixed subplot layout and consistent color mapping across conditions.
- Support one condition or all conditions.

Inputs:
- `results/integrated_samples/<Condition>/<Condition>_integrated.h5ad`

Outputs:
- `results/integrated_samples/panels/umap_panel_leiden.png`
- `results/integrated_samples/panels/umap_panel_SampleName.png`
- Panel figure includes embedded legend on the right side.
- Full color mapping CSV per key:
	- `results/integrated_samples/panels/umap_panel_<key>_legend.csv`

Run:
- List available condition folders:
	- `python scripts/RtoPython/plot_integrated_panels.py --list-conditions`
- Build panels for all conditions:
	- `python scripts/RtoPython/plot_integrated_panels.py --condition all`
- Build panels for one condition:
	- `python scripts/RtoPython/plot_integrated_panels.py --condition Normal`
- Custom color keys and layout:
	- `python scripts/RtoPython/plot_integrated_panels.py --condition all --color-keys "leiden,SampleName" --n-cols 4`
- Clean legend control:
	- `python scripts/RtoPython/plot_integrated_panels.py --condition all --legend-max-items 30 --legend-cols 2`

## Notes

- This conversion reproduces the QC visualization logic; it does not run full Seurat workflows.
- Additional R scripts (annotation/clustering/group-specific scripts) can be converted here incrementally.
