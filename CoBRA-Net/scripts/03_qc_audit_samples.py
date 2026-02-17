#!/usr/bin/env python3
"""03_qc_audit_samples.py: light QC audit over the selected cohort."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Run light QC audit on bigboss_chopped cohort.")
	parser.add_argument("--input", default="results/bigboss_chopped.csv", help="Path to chopped cohort CSV.")
	parser.add_argument("--data-dir", default="data/GSE161529_RAW", help="Path to raw matrix/barcode files.")
	parser.add_argument("--features", default="data/GSE161529_features.tsv", help="Path to features TSV.")
	parser.add_argument("--out-dir", default="results/qc_audit", help="Output folder for QC summary and plots.")
	parser.add_argument("--max-samples", type=int, default=0, help="Optional: limit number of samples for quick test.")
	parser.add_argument("--top-n-violin-samples", type=int, default=8, help="Number of largest samples to show in per-sample violin plots.")
	parser.add_argument("--max-cells-per-sample", type=int, default=5000, help="Max cells per sample to keep for violin plotting (for readability).")
	return parser.parse_args()


def main() -> int:
	args = parse_args()

	repo_root = Path(__file__).resolve().parents[1]
	if str(repo_root) not in sys.path:
		sys.path.insert(0, str(repo_root))

	from utils.qc_functions import (
		add_qc_flags,
		compute_basic_qc,
		extract_per_cell_qc,
		load_features,
		load_matrix_and_barcodes,
		make_cohort_plots,
		make_thesis_panel_figure,
		make_violin_plots,
	)

	in_path = (repo_root / args.input).resolve()
	data_dir = (repo_root / args.data_dir).resolve()
	features_path = (repo_root / args.features).resolve()
	out_dir = (repo_root / args.out_dir).resolve()
	out_dir.mkdir(parents=True, exist_ok=True)

	if not in_path.exists():
		print(f"ERROR: input cohort file not found: {in_path}")
		return 1

	cohort = pd.read_csv(in_path)
	if args.max_samples > 0:
		cohort = cohort.head(args.max_samples).copy()

	features = load_features(features_path)
	gene_names = features["gene_name"]

	rows = []
	per_cell_rows = []
	total = len(cohort)
	for idx, row in cohort.iterrows():
		sample_name = row["SampleName"]
		matrix_file = data_dir / row["MatrixFile"]
		barcodes_file = data_dir / row["BarcodesFile"]

		if not matrix_file.exists() or not barcodes_file.exists():
			print(f"[{idx + 1}/{total}] {sample_name}: missing files, skipping")
			continue

		matrix, barcodes = load_matrix_and_barcodes(matrix_file, barcodes_file)
		qc = compute_basic_qc(matrix, gene_names)
		per_cell = extract_per_cell_qc(matrix, gene_names)

		qc["barcode_count"] = int(len(barcodes))
		qc["barcode_match"] = int(matrix.shape[0] == len(barcodes))
		qc["SampleName"] = sample_name
		qc["GEO_ID"] = row.get("GEO_ID", "")
		qc["Condition"] = row.get("Condition", "")
		qc["Source"] = row.get("Source", "")
		qc["Gender"] = row.get("Gender", "")

		per_cell["SampleName"] = sample_name
		per_cell["Condition"] = row.get("Condition", "")
		if args.max_cells_per_sample > 0 and len(per_cell) > args.max_cells_per_sample:
			per_cell = per_cell.sample(n=args.max_cells_per_sample, random_state=42)
		per_cell_rows.append(per_cell)

		rows.append(qc)
		print(f"[{idx + 1}/{total}] {sample_name}: OK ({qc['n_cells']} cells, {qc['n_genes']} genes)")

	summary_df = pd.DataFrame(rows)
	if summary_df.empty:
		print("ERROR: no samples were processed.")
		return 1

	summary_df.to_csv(out_dir / "qc_audit_summary.csv", index=False)
	summary_with_flags = add_qc_flags(summary_df)
	summary_with_flags.to_csv(out_dir / "qc_audit_summary_with_flags.csv", index=False)
	make_cohort_plots(summary_df, out_dir)

	per_cell_df = pd.concat(per_cell_rows, ignore_index=True) if per_cell_rows else pd.DataFrame()
	if not per_cell_df.empty:
		make_violin_plots(
			per_cell_df=per_cell_df,
			summary_df=summary_df,
			out_dir=out_dir,
			top_n_samples=args.top_n_violin_samples,
		)
		make_thesis_panel_figure(per_cell_df=per_cell_df, out_dir=out_dir)

	print("\n=== QC audit complete ===")
	print(f"Samples processed: {len(summary_df)}")
	print(f"Summary: {out_dir / 'qc_audit_summary.csv'}")
	print(f"Flags:   {out_dir / 'qc_audit_summary_with_flags.csv'}")
	print(f"Plots:   {out_dir}")
	print("Violin plots: condition-level (all samples) + per-sample (largest subset)")
	print(f"Thesis panel: {out_dir / 'qc_thesis_panel_condition_violin.png'}")
	print("Note: this script audits quality and does not filter cells.")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
