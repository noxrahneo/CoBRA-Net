# Utils Overview

This folder contains reusable helper modules used by scripts in the CoBRA-Net pipeline.

## qc_functions.py

Purpose:
- Shared QC helpers used by `scripts/RtoPython/01_pre_qc_from_r.py` and
	`scripts/RtoPython/03_post_qc_from_r.py`.
- Loads features and matrix/barcode files.
- Computes basic per-sample QC metrics (counts, genes, sparsity, mito%).
- Adds simple QC flags and cohort-level plots.

Why it matters:
- Keeps analysis scripts lightweight and easy to read.
- Makes QC logic reusable and easier to maintain.

## Notes

- `data_loading.py` in the repo root is legacy and not used by the pipeline.
- `__pycache__/` appear automatically after running Python scripts.
- Utility files should contain reusable functions, not pipeline orchestration logic.
