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

1. `pre_qc_from_r.py`
2. `filtering_cells.py`
3. `post_qc_from_r.py`

See `scripts/RtoPython/README.md` for commands and outputs.
