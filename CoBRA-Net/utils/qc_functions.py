from __future__ import annotations

import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.io import mmread


def load_features(features_path: Path) -> pd.DataFrame:
    features = pd.read_csv(features_path, sep="\t", header=None)
    if features.shape[1] < 2:
        raise ValueError("Features file must have at least 2 columns.")
    features = features.iloc[:, :2].copy()
    features.columns = ["gene_id", "gene_name"]
    return features


def load_matrix_and_barcodes(matrix_path: Path, barcodes_path: Path):
    with gzip.open(matrix_path, "rb") as handle:
        matrix = mmread(handle).T.tocsr()  # cells x genes
    barcodes = pd.read_csv(barcodes_path, header=None)
    return matrix, barcodes


def compute_basic_qc(matrix, gene_names: pd.Series) -> dict:
    n_cells, n_genes = matrix.shape

    counts_per_cell = np.asarray(matrix.sum(axis=1)).ravel()
    genes_per_cell = np.asarray((matrix > 0).sum(axis=1)).ravel()
    cells_per_gene = np.asarray((matrix > 0).sum(axis=0)).ravel()

    mito_mask = gene_names.astype(str).str.upper().str.startswith("MT-").values
    if mito_mask.any():
        mito_counts = np.asarray(matrix[:, mito_mask].sum(axis=1)).ravel()
        pct_mito = np.divide(
            mito_counts,
            counts_per_cell,
            out=np.zeros_like(mito_counts, dtype=float),
            where=counts_per_cell > 0,
        ) * 100.0
    else:
        pct_mito = np.full(n_cells, np.nan)

    sparsity = 1.0 - (matrix.nnz / (n_cells * n_genes))

    return {
        "n_cells": int(n_cells),
        "n_genes": int(n_genes),
        "total_counts": float(counts_per_cell.sum()),
        "mean_counts_per_cell": float(np.mean(counts_per_cell)),
        "median_counts_per_cell": float(np.median(counts_per_cell)),
        "mean_genes_per_cell": float(np.mean(genes_per_cell)),
        "median_genes_per_cell": float(np.median(genes_per_cell)),
        "pct_mito_mean": float(np.nanmean(pct_mito)) if not np.isnan(pct_mito).all() else np.nan,
        "pct_mito_median": float(np.nanmedian(pct_mito)) if not np.isnan(pct_mito).all() else np.nan,
        "genes_detected_in_ge_3_cells": int((cells_per_gene >= 3).sum()),
        "sparsity_pct": float(sparsity * 100.0),
    }


def extract_per_cell_qc(matrix, gene_names: pd.Series) -> pd.DataFrame:
    counts_per_cell = np.asarray(matrix.sum(axis=1)).ravel()
    genes_per_cell = np.asarray((matrix > 0).sum(axis=1)).ravel()

    mito_mask = gene_names.astype(str).str.upper().str.startswith("MT-").values
    if mito_mask.any():
        mito_counts = np.asarray(matrix[:, mito_mask].sum(axis=1)).ravel()
        pct_mito = np.divide(
            mito_counts,
            counts_per_cell,
            out=np.zeros_like(mito_counts, dtype=float),
            where=counts_per_cell > 0,
        ) * 100.0
    else:
        pct_mito = np.full(matrix.shape[0], np.nan)

    return pd.DataFrame(
        {
            "n_counts": counts_per_cell,
            "n_genes": genes_per_cell,
            "pct_mito": pct_mito,
        }
    )


def add_qc_flags(summary_df: pd.DataFrame) -> pd.DataFrame:
    flags = []
    for _, row in summary_df.iterrows():
        issues = []
        if row["median_genes_per_cell"] < 300:
            issues.append("low_genes_per_cell")
        if pd.notna(row["pct_mito_median"]) and row["pct_mito_median"] > 20:
            issues.append("high_mito")
        if row.get("barcode_match", 1) == 0:
            issues.append("barcode_mismatch")
        flags.append(",".join(issues) if issues else "ok")

    out = summary_df.copy()
    out["qc_flag"] = flags
    return out


def make_cohort_plots(summary_df: pd.DataFrame, out_dir: Path) -> None:
    sns.set(style="whitegrid")

    plt.figure(figsize=(12, 5))
    sns.barplot(data=summary_df.sort_values("n_cells", ascending=False), x="SampleName", y="n_cells")
    plt.xticks(rotation=90)
    plt.title("Cells per sample (QC audit)")
    plt.tight_layout()
    plt.savefig(out_dir / "qc_cells_per_sample.png", dpi=150)
    plt.close()

    plt.figure(figsize=(12, 5))
    sns.barplot(
        data=summary_df.sort_values("median_genes_per_cell", ascending=False),
        x="SampleName",
        y="median_genes_per_cell",
    )
    plt.xticks(rotation=90)
    plt.title("Median genes per cell (QC audit)")
    plt.tight_layout()
    plt.savefig(out_dir / "qc_median_genes_per_cell.png", dpi=150)
    plt.close()

    if "Condition" in summary_df.columns:
        plt.figure(figsize=(10, 5))
        sns.boxplot(data=summary_df, x="Condition", y="pct_mito_median")
        plt.xticks(rotation=45, ha="right")
        plt.title("Median mitochondrial % by condition (QC audit)")
        plt.tight_layout()
        plt.savefig(out_dir / "qc_mito_by_condition.png", dpi=150)
        plt.close()


def make_violin_plots(
    per_cell_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    out_dir: Path,
    top_n_samples: int = 8,
) -> None:
    if per_cell_df.empty:
        return

    # Full cohort: violins by condition
    if "Condition" in per_cell_df.columns:
        metrics = ["n_counts", "n_genes", "pct_mito"]
        for metric in metrics:
            plot_df = per_cell_df.copy()
            if metric == "pct_mito":
                plot_df = plot_df[plot_df[metric].notna()].copy()
                if plot_df.empty:
                    continue

            plt.figure(figsize=(12, 6))
            sns.violinplot(data=plot_df, x="Condition", y=metric, inner="quartile", cut=0)
            plt.xticks(rotation=45, ha="right")
            plt.title(f"{metric} distribution by condition")
            plt.tight_layout()
            plt.savefig(out_dir / f"qc_violin_condition_{metric}.png", dpi=150)
            plt.close()

    # Subset: per-sample violins for largest samples only
    if "SampleName" in per_cell_df.columns and "n_cells" in summary_df.columns:
        top_samples = (
            summary_df.sort_values("n_cells", ascending=False)["SampleName"]
            .head(top_n_samples)
            .tolist()
        )
        sample_df = per_cell_df[per_cell_df["SampleName"].isin(top_samples)].copy()

        if not sample_df.empty:
            metrics = ["n_counts", "n_genes", "pct_mito"]
            for metric in metrics:
                plot_df = sample_df.copy()
                if metric == "pct_mito":
                    plot_df = plot_df[plot_df[metric].notna()].copy()
                    if plot_df.empty:
                        continue

                plt.figure(figsize=(12, 6))
                sns.violinplot(data=plot_df, x="SampleName", y=metric, inner="quartile", cut=0)
                plt.xticks(rotation=45, ha="right")
                plt.title(f"{metric} distribution for top {len(top_samples)} largest samples")
                plt.tight_layout()
                plt.savefig(out_dir / f"qc_violin_top_samples_{metric}.png", dpi=150)
                plt.close()


def make_thesis_panel_figure(per_cell_df: pd.DataFrame, out_dir: Path) -> None:
    if per_cell_df.empty or "Condition" not in per_cell_df.columns:
        return

    panel_metrics = ["n_counts", "n_genes", "pct_mito"]
    available_metrics = [metric for metric in panel_metrics if metric in per_cell_df.columns]
    if not available_metrics:
        return

    fig, axes = plt.subplots(1, len(available_metrics), figsize=(6 * len(available_metrics), 5))
    if len(available_metrics) == 1:
        axes = [axes]

    for ax, metric in zip(axes, available_metrics):
        plot_df = per_cell_df.copy()
        if metric == "pct_mito":
            plot_df = plot_df[plot_df[metric].notna()].copy()
            if plot_df.empty:
                ax.set_visible(False)
                continue

        sns.violinplot(data=plot_df, x="Condition", y=metric, inner="quartile", cut=0, ax=ax)
        ax.set_title(f"{metric} by condition")
        ax.tick_params(axis="x", rotation=45)

    fig.suptitle("QC audit panel (cohort-level)", y=1.02)
    fig.tight_layout()
    fig.savefig(out_dir / "qc_thesis_panel_condition_violin.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
