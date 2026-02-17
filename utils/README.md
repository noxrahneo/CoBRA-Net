# Utils Overview

This folder contains reusable helper modules used by scripts in the CoBRA-Net pipeline.

## data_loading.py

Purpose:
- Shared functions for loading metadata and sample files.
- Helpers for locating matrix/barcode files from inventory tables.
- Common data-loading logic used by downstream scripts.

Why it matters:
- Keeps file I/O logic in one place.
- Reduces duplication across scripts.

## qc_functions.py

Purpose:
- Shared QC-audit helpers used by `scripts/03_qc_audit_samples.py`.
- Loads features and matrix/barcode files.
- Computes basic per-sample QC metrics (counts, genes, sparsity, mito%).
- Adds simple QC flags and cohort-level plots.

Why it matters:
- Keeps analysis scripts lightweight and easy to read.
- Makes QC logic reusable and easier to maintain.

## Notes

- `__pycache__/` appear automatically after running Python scripts.
- Utility files should contain reusable functions, not pipeline orchestration logic.
