# R to Python Conversions

This folder contains Python conversions of selected scripts from `HumanBreast10X-main/RCode`.

## Why this folder exists

- Keep direct R→Python translations in one place.
- Preserve the logic of the original publication code.
- Make conversions easy to review, compare and cite in thesis methods.

## Current conversions

Note:
- Non-translation analysis utilities moved to `scripts/analysis/`.
- R→Python translations remain in this folder only.
- Translation marker: scripts include `R_TRANSLATED: yes` in module docstrings.

Translated scripts in execution order:
- `01_pre_qc_from_r.py`
- `02_filtering_cells.py`
- `03_post_qc_from_r.py`
- `04_per_sample_preprocess.py`
- `05_validate_preprocess.py`
- `06_integrate_per_condition.py`
- `07_plot_integrated_results.py`
- `08_plot_integrated_panels.py`
- `09_annotate_clusters.py`

## Quick run (all conditions)

From repository root, run the post-annotation analyses in this order:

```bash
python3 -m pip install --user gseapy statsmodels
python3 scripts/analysis/04_marker_heatmaps.py --condition all --group-col cell_type_annot
python3 scripts/analysis/05_composition_tests.py --condition all --group-col cell_type_annot --sample-col SampleName
python3 scripts/analysis/06_kegg_enrichment.py --condition all
```

Optional cluster-level variant:

```bash
python3 scripts/analysis/04_marker_heatmaps.py --condition all --group-col leiden
python3 scripts/analysis/05_composition_tests.py --condition all --group-col leiden --sample-col SampleName
```

## Stage warehouse logs

Each stage now writes a `warehouse.csv` file (append-only) to the stage output
directory for run provenance.

Examples:
- Annotation: `results/stages/04_annotation/warehouse.csv`
- Composition: `results/stages/05_composition/warehouse.csv`
- KEGG: `results/stages/06_kegg/warehouse.csv`
- Thesis pack: `results/stages/90_reports/thesis_pack/warehouse.csv`

Warehouse columns:
- `input_file`: stage input file/path.
- `output_file`: primary produced file.
- `script`: script path used.
- `date_utc`: UTC timestamp of the run record.
- `params_hash`: short hash of CLI parameter set.
- `condition`: condition label (or `all`).
- `stage`: pipeline stage name.

### `01_pre_qc_from_r.py`

Purpose:
- Runs pre-filter QC directly on raw matrices from `results/bigboss_chopped.csv`.
- Computes sample-level and per-cell QC metrics before filtering.

Outputs:
- `results/stages/01_qc/qc_pre/pre_qc_summary.csv`
- `results/stages/01_qc/qc_pre/pre_qc_per_cell.parquet` (or CSV fallback)
- `results/stages/01_qc/qc_pre/pre_qc_boxplots.png`
- `results/stages/01_qc/qc_pre/pre_qc_violins.png`

Run:
- `python scripts/RtoPython/01_pre_qc_from_r.py`
- `python scripts/RtoPython/01_pre_qc_from_r.py --condition "ER+ tumor"`

### `02_filtering_cells.py`

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
- `results/stages/02_preprocess/filtered/<Condition>/<SampleName>_filtered.h5ad`
- `results/stages/02_preprocess/filtered/filtering_summary.csv`

Run:
- `python scripts/RtoPython/02_filtering_cells.py`
- `python scripts/RtoPython/02_filtering_cells.py --condition Normal`
- `python scripts/RtoPython/02_filtering_cells.py --list-conditions`

Important clarification:
- This script performs **actual filtering** and writes new filtered outputs.
- It is still a scoped conversion focused on filtering logic, not a full reproduction of every downstream Seurat step in `NormTotal.R`.

### `03_post_qc_from_r.py`

Purpose:
- Runs post-filter QC on filtered `.h5ad` files from `results/stages/02_preprocess/filtered`.
- Uses the same metric definitions as pre-QC.
- Compares pre/post QC when pre summary is available.

Outputs:
- `results/stages/01_qc/qc_post/post_qc_summary.csv`
- `results/stages/01_qc/qc_post/post_qc_per_cell.parquet` (or CSV fallback)
- `results/stages/01_qc/qc_post/post_qc_boxplots.png`
- `results/stages/01_qc/qc_post/post_qc_violins.png`
- `results/stages/01_qc/qc_post/pre_post_qc_comparison.csv` (if pre summary exists)
- `results/stages/01_qc/qc_post/pre_post_qc_scatter.png` (if pre summary exists)
- `results/stages/01_qc/qc_post/pre_post_qc_violins.png` (if pre per-cell exists)

Run:
- `python scripts/RtoPython/03_post_qc_from_r.py`
- `python scripts/RtoPython/03_post_qc_from_r.py --condition "ER+ tumor"`

### `04_per_sample_preprocess.py`

Purpose:
- Run per-sample preprocessing before integration.
- Seurat-style step mapping: NormalizeData -> FindVariableFeatures -> ScaleData.
- Support one condition or all conditions from condition folders.

Outputs:
- Preprocessed sample files:
	- `results/stages/02_preprocess/preprocessed/<Condition>/<SampleName>_preprocessed.h5ad`
- Summary CSV:
	- one condition: `results/stages/02_preprocess/preprocessed/<Condition>/preprocess_summary.csv`
	- all conditions: `results/stages/02_preprocess/preprocessed/preprocess_summary.csv`

Run:
- List available condition folders:
	- `python scripts/RtoPython/04_per_sample_preprocess.py --list-conditions`
- Run one condition:
	- `python scripts/RtoPython/04_per_sample_preprocess.py --condition Normal`
- Run one condition by label (safe folder mapping supported):
	- `python scripts/RtoPython/04_per_sample_preprocess.py --condition "ER+ tumor"`
- Run all conditions:
	- `python scripts/RtoPython/04_per_sample_preprocess.py --condition all`

### `05_validate_preprocess.py`

Purpose:
- Validate per-sample preprocessing outputs from `04_per_sample_preprocess.py`.
- Compare filtered vs preprocessed files sample-by-sample.
- Check preprocessing artifacts and consistency:
	- `counts` layer exists
	- `highly_variable` column exists
	- HVG count passes selected mode (`strict` or `relaxed`)
	- no non-finite values in expression matrix
	- cell/gene dimensions are preserved

Outputs:
- Per-condition CSV:
	- `results/stages/02_preprocess/preprocess_validation/preprocess_validation_<Condition>.csv`
- Run-level summary JSON:
	- `results/stages/02_preprocess/preprocess_validation/preprocess_validation_<Condition>.json`
- If `--condition all`, also writes:
	- `results/stages/02_preprocess/preprocess_validation/preprocess_validation_all.csv`
	- `results/stages/02_preprocess/preprocess_validation/preprocess_validation_all.json`

Run:
- List available condition folders:
	- `python scripts/RtoPython/05_validate_preprocess.py --list-conditions`
- Validate one condition:
	- `python scripts/RtoPython/05_validate_preprocess.py --condition Normal`
- Validate one condition (strict HVG mode):
	- `python scripts/RtoPython/05_validate_preprocess.py --condition Normal --hvg-mode strict --expected-hvg 3000`
- Validate one condition (relaxed HVG mode):
	- `python scripts/RtoPython/05_validate_preprocess.py --condition Normal --hvg-mode relaxed --hvg-tolerance 0.05 --expected-hvg 3000`
- Validate all conditions:
	- `python scripts/RtoPython/05_validate_preprocess.py --condition all --hvg-mode relaxed`

### `06_integrate_per_condition.py`

Purpose:
- Integrate preprocessed samples per condition into a joint condition object.
- Seurat-style mapping:
	- FindIntegrationAnchors / IntegrateData (Python equivalent)
	- joint PCA / neighbors / UMAP / Leiden
- Supports one condition or all conditions.

Outputs:
- Integrated condition object:
	- `results/stages/03_integration/integrated/<Condition>/<Condition>_integrated.h5ad`
- Per-cell metadata table:
	- `results/stages/03_integration/integrated/<Condition>/<Condition>_cell_metadata.csv`
- UMAP coordinates table:
	- `results/stages/03_integration/integrated/<Condition>/<Condition>_umap.csv`
- Integration summary:
	- one condition: `results/stages/03_integration/integration_summary_<Condition>.csv`
	- all conditions: `results/stages/03_integration/integration_summary.csv`

Run:
- List available condition folders:
	- `python scripts/RtoPython/06_integrate_per_condition.py --list-conditions`
- Integrate one condition:
	- `python scripts/RtoPython/06_integrate_per_condition.py --condition Normal`
- Integrate one condition by label (safe folder mapping supported):
	- `python scripts/RtoPython/06_integrate_per_condition.py --condition "ER+ tumor"`
- Integrate all conditions:
	- `python scripts/RtoPython/06_integrate_per_condition.py --condition all`

### `07_plot_integrated_results.py`

Purpose:
- Create post-integration visualization plots per condition.
- Save UMAP plots by cluster (`leiden`) and by sample (`SampleName`).
- Optional t-SNE plotting for article-style embedding comparison.

Inputs:
- `results/stages/03_integration/integrated/<Condition>/<Condition>_integrated.h5ad`

Outputs:
- `results/stages/03_integration/integrated/<Condition>/figures/<Condition>_umap_leiden.png`
- `results/stages/03_integration/integrated/<Condition>/figures/<Condition>_umap_sample.png`
- Optional (`--compute-tsne`):
	- `results/stages/03_integration/integrated/<Condition>/figures/<Condition>_tsne_leiden.png`
	- `results/stages/03_integration/integrated/<Condition>/figures/<Condition>_tsne_sample.png`

Run:
- List available condition folders:
	- `python scripts/RtoPython/07_plot_integrated_results.py --list-conditions`
- Plot one condition:
	- `python scripts/RtoPython/07_plot_integrated_results.py --condition Normal`
- Plot all conditions:
	- `python scripts/RtoPython/07_plot_integrated_results.py --condition all`
- Plot all conditions with t-SNE:
	- `python scripts/RtoPython/07_plot_integrated_results.py --condition all --compute-tsne`
- If legends are crowded/cut, set legend options:
	- `python scripts/RtoPython/07_plot_integrated_results.py --condition Normal --legend-loc "right margin" --legend-fontsize 7`

### `08_plot_integrated_panels.py`

Purpose:
- Build publication-style multi-condition UMAP panel figures.
- Use fixed subplot layout and consistent color mapping across conditions.
- Support one condition or all conditions.

Inputs:
- `results/stages/03_integration/integrated/<Condition>/<Condition>_integrated.h5ad`

Outputs:
- `results/stages/03_integration/integrated/panels/umap_panel_leiden.png`
- `results/stages/03_integration/integrated/panels/umap_panel_SampleName.png`
- Panel figure includes embedded legend on the right side.
- Full color mapping CSV per key:
	- `results/stages/03_integration/integrated/panels/umap_panel_<key>_legend.csv`

Run:
- List available condition folders:
	- `python scripts/RtoPython/08_plot_integrated_panels.py --list-conditions`
- Build panels for all conditions:
	- `python scripts/RtoPython/08_plot_integrated_panels.py --condition all`
- Build panels for one condition:
	- `python scripts/RtoPython/08_plot_integrated_panels.py --condition Normal`
- Custom color keys and layout:
	- `python scripts/RtoPython/08_plot_integrated_panels.py --condition all --color-keys "leiden,SampleName" --n-cols 4`
- Clean legend control:
	- `python scripts/RtoPython/08_plot_integrated_panels.py --condition all --legend-cols 3`
	- (optional limit) `python scripts/RtoPython/08_plot_integrated_panels.py --condition all --legend-max-items 30 --legend-cols 3`

### `09_annotate_clusters.py`

Purpose:
- Annotate clusters for one integrated condition or all integrated conditions.
- Build marker/DE support per cluster.
- Score signatures and infer cluster -> cell-type labels.
- Export mapping tables and regenerated UMAPs colored by annotation.

Inputs:
- `results/stages/03_integration/integrated/<Condition>/<Condition>_integrated.h5ad`
- Optional signature file (default):
	- `data/HumanBreast10X-main/Signatures/ImmuneMarkers2.txt`

Outputs:
- `results/stages/04_annotation/<Condition>/<Condition>_annotated.h5ad`
- `results/stages/04_annotation/<Condition>/<Condition>_cluster_markers_top.csv`
- `results/stages/04_annotation/<Condition>/<Condition>_cluster_signature_scores.csv`
- `results/stages/04_annotation/<Condition>/<Condition>_cluster_to_celltype_mapping.csv`
- `results/stages/04_annotation/<Condition>/<Condition>_signature_gene_coverage.csv`
- Optional confusion report:
	- `results/stages/04_annotation/<Condition>/<Condition>_confusion_previous_vs_annot.csv`
- Annotation summary:
	- one condition: `results/stages/04_annotation/<Condition>/<Condition>_annotation_summary.csv`
	- all conditions: `results/stages/04_annotation/annotation_summary.csv`
- Annotated UMAP figures:
	- `results/stages/04_annotation/<Condition>/figures/<Condition>_umap_cell_type_annot.png`
	- `results/stages/04_annotation/<Condition>/figures/<Condition>_umap_leiden.png`

Run:
- List available condition folders:
	- `python scripts/RtoPython/09_annotate_clusters.py --list-conditions`
- Normal-first run:
	- `python scripts/RtoPython/09_annotate_clusters.py --condition Normal`
- Run all conditions:
	- `python scripts/RtoPython/09_annotate_clusters.py --condition all`
- With previous labels for confusion report:
	- `python scripts/RtoPython/09_annotate_clusters.py --condition Normal --previous-label-col cell_type`

### `04_marker_heatmaps.py`

Purpose:
- Build marker heatmaps from annotated outputs.
- Use per-cluster marker table + annotated object to compute group-level mean expression.
- Export both expression matrix and row z-score matrix used for plotting.

Inputs:
- `results/stages/04_annotation/<Condition>/<Condition>_annotated.h5ad`
- `results/stages/04_annotation/<Condition>/<Condition>_cluster_markers_top.csv`

Outputs:
- `results/stages/04_annotation/<Condition>/<Condition>_marker_heatmap_expression.csv`
- `results/stages/04_annotation/<Condition>/<Condition>_marker_heatmap_zscore.csv`
- `results/stages/04_annotation/<Condition>/figures/<Condition>_marker_heatmap_<group_col>.png`
- If `--condition all`:
	- `results/stages/04_annotation/marker_heatmap_summary.csv`

Run:
- List available annotation condition folders:
	- `python scripts/analysis/04_marker_heatmaps.py --list-conditions`
- Run all conditions (default):
	- `python scripts/analysis/04_marker_heatmaps.py --condition all`
- Run one condition:
	- `python scripts/analysis/04_marker_heatmaps.py --condition Normal`
- Build heatmap by clusters instead of cell-type labels:
	- `python scripts/analysis/04_marker_heatmaps.py --condition all --group-col leiden`

### `05_composition_tests.py`

Purpose:
- Test global cell composition shifts across conditions.
- Build sample-level count table by cell group (`cell_type_annot` by default).
- Fit quasi-Poisson style GLM (Poisson with overdispersion scaling).

Inputs:
- `results/stages/04_annotation/<Condition>/<Condition>_annotated.h5ad`

Outputs:
- `results/stages/05_composition/composition_counts_long.csv`
- `results/stages/05_composition/composition_proportions_by_sample.csv`
- `results/stages/05_composition/composition_glm_global_test.csv`
- `results/stages/05_composition/composition_glm_coefficients.csv`
- `results/stages/05_composition/composition_condition_stacked_bar.png`

Run:
- List available annotation condition folders:
	- `python scripts/analysis/05_composition_tests.py --list-conditions`
- Run global test across all conditions:
	- `python scripts/analysis/05_composition_tests.py --condition all`
- Test using clusters instead of annotated labels:
	- `python scripts/analysis/05_composition_tests.py --condition all --group-col leiden`

Dependency:
- Requires `statsmodels` for GLM fitting.

### `06_kegg_enrichment.py`

Purpose:
- Run KEGG enrichment from cluster marker genes.
- Process each condition and cluster independently.
- Save enrichment tables and skipped-cluster logs.

Inputs:
- `results/stages/04_annotation/<Condition>/<Condition>_cluster_markers_top.csv`

Outputs:
- Per condition:
	- `results/stages/06_kegg/<Condition>/<Condition>_kegg_enrichment.csv`
	- `results/stages/06_kegg/<Condition>/<Condition>_kegg_skipped.csv`
	- `results/stages/06_kegg/<Condition>/<Condition>_kegg_interpretation.csv`
- Combined (when available):
	- `results/stages/06_kegg/kegg_enrichment_all.csv`
	- `results/stages/06_kegg/kegg_skipped_all.csv`
	- `results/stages/06_kegg/kegg_interpretation_all.csv`

Run:
- List available annotation condition folders:
	- `python scripts/analysis/06_kegg_enrichment.py --list-conditions`
- Run all conditions:
	- `python scripts/analysis/06_kegg_enrichment.py --condition all`
- Run one condition:
	- `python scripts/analysis/06_kegg_enrichment.py --condition Normal`
- Use alternative Enrichr library:
	- `python scripts/analysis/06_kegg_enrichment.py --condition all --gene-set KEGG_2019_Human`
- Tune interpretation filter for thesis-safe reporting:
	- `python scripts/analysis/06_kegg_enrichment.py --condition all --interpret-padj 0.05 --interpret-min-overlap 3 --interpret-max-terms 15`

Dependency:
- Requires `gseapy` for Enrichr KEGG queries.

Interpretation note:
- `*_kegg_enrichment.csv` is the raw enrichment table.
- `*_kegg_interpretation.csv` is a filtered, thesis-oriented subset:
	- keeps terms with `Adjusted P-value <= --interpret-padj`
	- keeps terms with overlap hits `>= --interpret-min-overlap`
	- keeps top terms per cluster (`--interpret-max-terms`)
- Disease terms indicate pathway-level similarity of marker genes, not proof of causal disease state in your sample.

## Output column reference (CSV headers)

This section explains the main CSV headers produced by the current pipeline.

### Preprocess summary

File:
- `results/stages/02_preprocess/preprocessed/<Condition>/preprocess_summary.csv`

Headers:
- `Condition`: condition folder name.
- `SampleName`: sample identifier.
- `InputFile`: source filtered `.h5ad` file used for preprocessing.
- `OutputFile`: generated preprocessed `.h5ad` file.
- `Cells`: number of cells in the sample.
- `Genes`: number of genes before HVG selection.
- `HVGs`: number of genes flagged as highly variable.
- `TargetSum`: normalization target used in `normalize_total`.
- `NTopGenes`: requested HVG count.
- `ScaleMaxValue`: clipping value used in scaling.

### Preprocess validation report

File:
- `results/stages/02_preprocess/preprocess_validation/preprocess_validation_<Condition>.csv`

Headers:
- `SampleName`: sample identifier.
- `filtered_file`: filtered input `.h5ad`.
- `preprocessed_file`: expected preprocessed `.h5ad`.
- `exists_preprocessed`: whether preprocessed file exists.
- `expected_hvg`: requested HVG target for validation.
- `hvg_mode`: `strict` or `relaxed` validation mode.
- `hvg_tolerance`: tolerance fraction for `relaxed` mode.
- `cells_filtered`, `cells_preprocessed`: cell counts before/after preprocessing.
- `genes_filtered`, `genes_preprocessed`: gene counts before/after preprocessing.
- `has_counts_layer`: whether `counts` layer exists.
- `has_hvg_column`: whether `highly_variable` column exists.
- `hvg_count`: observed HVG count.
- `hvg_target`: effective target (min of expected HVG and n_vars).
- `x_is_finite`: whether matrix values are finite (no NaN/Inf).
- `status`: `PASS` or `FAIL`.
- `issues`: semicolon-separated validation issues.
- `Condition`: condition name.

### Integration summary

File:
- `results/stages/03_integration/integration_summary_<Condition>.csv`

Headers:
- `Condition`: integrated condition.
- `Samples`: number of samples merged for that condition.
- `Cells`: number of cells in integrated object.
- `GenesAfterHVG`: genes retained after HVG filtering.
- `Clusters`: number of Leiden clusters.
- `IntegrationMethod`: integration method used (`combat` or `none`).
- `IntegratedFile`: integrated `.h5ad` path.
- `CellMetadata`: per-cell metadata CSV path.
- `UmapFile`: UMAP coordinate CSV path.

### Integrated per-cell metadata

File:
- `results/stages/03_integration/integrated/<Condition>/<Condition>_cell_metadata.csv`

Headers:
- `Unnamed: 0`: cell barcode/cell ID (CSV index).
- `SampleName`: sample identifier for each cell.
- `Condition`: condition label for each cell.
- `leiden`: cluster assignment.

### Integrated UMAP coordinates

File:
- `results/stages/03_integration/integrated/<Condition>/<Condition>_umap.csv`

Headers:
- `Unnamed: 0`: cell barcode/cell ID (CSV index).
- `UMAP1`, `UMAP2`: UMAP embedding coordinates.

### Panel legend mapping

File:
- `results/stages/03_integration/integrated/panels/umap_panel_<key>_legend.csv`

Headers:
- `category`: category label shown in panel coloring.
- `color_hex`: hex color code assigned to that category.

### Annotation markers

File:
- `results/stages/04_annotation/<Condition>/<Condition>_cluster_markers_top.csv`

Headers:
- `cluster`: cluster ID.
- `names`: gene name.
- `scores`: rank score from Scanpy differential test.
- `logfoldchanges`: log fold-change for cluster vs reference.
- `pvals`: raw p-value.
- `pvals_adj`: multiple-testing adjusted p-value.

### Annotation signature scores (cluster-level)

File:
- `results/stages/04_annotation/<Condition>/<Condition>_cluster_signature_scores.csv`

Headers:
- `leiden`: cluster ID.
- `sig_*`: mean signature score for each signature in that cluster.

### Cluster-to-cell-type mapping

File:
- `results/stages/04_annotation/<Condition>/<Condition>_cluster_to_celltype_mapping.csv`

Headers:
- `leiden`: cluster ID.
- `predicted_cell_type`: inferred label from highest signature score.
- `score_top1`, `score_top2`: top two signature scores.
- `score_margin`: `score_top1 - score_top2` (confidence proxy).
- `top1_signature`, `top2_signature`: top scoring signatures.
- `marker_overlap_n`: overlap count between top markers and predicted signature genes.
- `marker_overlap_genes`: overlapping genes (`;`-separated).

### Signature coverage

File:
- `results/stages/04_annotation/<Condition>/<Condition>_signature_gene_coverage.csv`

Headers:
- `signature`: signature name.
- `genes_present_n`: number of signature genes found in the dataset.

### Annotation summary

File:
- `results/stages/04_annotation/<Condition>/<Condition>_annotation_summary.csv`

Headers:
- `condition`: condition name.
- `cells`: total cells in annotated object.
- `genes`: total genes in annotated object.
- `clusters`: number of clusters used for annotation.
- `annotated_cell_types`: number of inferred cell-type labels.
- `integrated_input`: source integrated `.h5ad` path.
- `annotated_output`: annotated `.h5ad` path.

### KEGG interpretation table

File:
- `results/stages/06_kegg/<Condition>/<Condition>_kegg_interpretation.csv`

Headers:
- `condition`: condition name.
- `cluster`: cluster ID used for marker set.
- `Term`: KEGG pathway/disease term name.
- `Adjusted P-value`: multiple-testing adjusted enrichment p-value (FDR-like).
- `P-value`: raw enrichment p-value.
- `Overlap`: overlap string (`hits/term_size`).
- `overlap_hits`: overlapping marker-gene count.
- `term_gene_count`: total genes in the KEGG term.
- `overlap_ratio`: `overlap_hits / term_gene_count`.
- `Genes`: overlapping genes from your marker list.
- `input_genes_n`: number of marker genes submitted for that cluster.
- `thesis_interpretation`: thesis-safe interpretation sentence.

Biological meaning:
- Enrichment indicates over-represented programs among marker genes.
- Disease terms represent pathway-level similarity, not direct proof of causality in your sample.

## Notes

- This conversion reproduces the QC visualization logic; it does not run full Seurat workflows.
- Additional R scripts (annotation/clustering/group-specific scripts) can be converted here incrementally.
