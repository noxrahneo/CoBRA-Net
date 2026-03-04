#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot KEGG enrichment results from merged CSV tables"
    )
    parser.add_argument(
        "--input-file",
        default="results/stages/06_kegg/kegg_interpretation_all.csv",
        help="Merged KEGG interpretation CSV",
    )
    parser.add_argument(
        "--output-dir",
        default="results/stages/06_kegg/figures",
        help="Directory for KEGG figure outputs",
    )
    parser.add_argument(
        "--top-terms-per-condition",
        type=int,
        default=15,
        help="Number of top rows (lowest adjusted p-value) per condition",
    )
    parser.add_argument(
        "--top-terms-heatmap",
        type=int,
        default=25,
        help="Number of globally top terms to include in condition heatmap",
    )
    parser.add_argument("--dpi", type=int, default=180, help="Figure DPI")
    return parser.parse_args()


def resolve_base(path_like: str) -> Path:
    path = Path(path_like)
    return path if path.is_absolute() else REPO_ROOT / path


def clean_condition_name(name: str) -> str:
    return (
        name.replace(" ", "_")
        .replace("/", "-")
        .replace("__", "_")
        .replace("___", "_")
    )


def ensure_numeric(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["Adjusted P-value"] = pd.to_numeric(
        out["Adjusted P-value"], errors="coerce"
    )
    out["overlap_hits"] = pd.to_numeric(out["overlap_hits"], errors="coerce")
    out["cluster"] = out["cluster"].astype(str)
    out = out.dropna(subset=["Adjusted P-value"])
    out = out[out["Adjusted P-value"] > 0]
    out["neglog10_fdr"] = -np.log10(out["Adjusted P-value"].clip(lower=1e-300))
    return out


def plot_condition_dotplot(
    condition_df: pd.DataFrame,
    condition: str,
    out_file: Path,
    top_n: int,
    dpi: int,
) -> None:
    top = condition_df.nsmallest(top_n, "Adjusted P-value").copy()
    top = top.sort_values("Adjusted P-value", ascending=False)
    top["label"] = top["Term"].astype(str) + " (C" + top["cluster"] + ")"

    fig_h = max(5.5, 0.42 * len(top))
    fig, ax = plt.subplots(figsize=(10.5, fig_h))

    cluster_codes = pd.Categorical(top["cluster"]).codes
    scatter = ax.scatter(
        top["neglog10_fdr"],
        top["label"],
        s=(top["overlap_hits"].fillna(1).clip(lower=1) * 24),
        c=cluster_codes,
        cmap="tab20",
        alpha=0.85,
        edgecolors="black",
        linewidths=0.25,
    )

    ax.set_xlabel("-log10(Adjusted P-value)")
    ax.set_ylabel("KEGG term (cluster)")
    ax.set_title(f"KEGG enrichment top terms: {condition}")
    ax.grid(axis="x", linestyle="--", alpha=0.35)

    handles, labels = scatter.legend_elements(prop="colors", alpha=0.85)
    cluster_labels = list(pd.Categorical(top["cluster"]).categories)
    if handles and cluster_labels:
        ax.legend(
            handles,
            cluster_labels,
            title="Cluster",
            loc="lower right",
            fontsize=8,
            title_fontsize=9,
            frameon=True,
        )

    plt.tight_layout()
    fig.savefig(out_file, dpi=dpi, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)


def plot_condition_term_heatmap(
    df: pd.DataFrame,
    out_file: Path,
    top_terms: int,
    dpi: int,
) -> None:
    best_by_term = (
        df.groupby("Term", observed=False)["Adjusted P-value"]
        .min()
        .sort_values()
    )
    keep_terms = best_by_term.head(top_terms).index

    hdf = df[df["Term"].isin(keep_terms)].copy()
    hdf = hdf.groupby(
        ["condition", "Term"], observed=False
    )["neglog10_fdr"].max()
    heat = hdf.unstack("Term").fillna(0.0)

    if heat.empty:
        return

    conditions = list(heat.index)
    terms = list(heat.columns)

    fig_w = max(10, 0.45 * len(terms))
    fig_h = max(4.5, 0.45 * len(conditions))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    img = ax.imshow(heat.values, aspect="auto", cmap="viridis")

    ax.set_xticks(range(len(terms)))
    ax.set_xticklabels(terms, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(conditions)))
    ax.set_yticklabels(conditions, fontsize=9)
    ax.set_title("KEGG enrichment heatmap by condition")

    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label("max -log10(Adjusted P-value)")

    plt.tight_layout()
    fig.savefig(out_file, dpi=dpi, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)


def main() -> None:
    args = parse_args()

    input_file = resolve_base(args.input_file)
    out_dir = resolve_base(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not input_file.exists():
        raise FileNotFoundError(
            f"Missing KEGG interpretation file: {input_file}"
        )

    df = pd.read_csv(input_file)
    required = {
        "condition",
        "cluster",
        "Term",
        "Adjusted P-value",
        "overlap_hits",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            "Input file is missing required columns: "
            + ", ".join(sorted(missing))
        )

    df = ensure_numeric(df)
    if df.empty:
        raise ValueError("No usable KEGG rows found after numeric filtering")

    for condition, cdf in df.groupby("condition", observed=False):
        stem = clean_condition_name(str(condition))
        out_file = out_dir / f"kegg_{stem}_top_terms_dotplot.png"
        plot_condition_dotplot(
            condition_df=cdf,
            condition=str(condition),
            out_file=out_file,
            top_n=args.top_terms_per_condition,
            dpi=args.dpi,
        )

    heatmap_file = out_dir / "kegg_conditions_top_terms_heatmap.png"
    plot_condition_term_heatmap(
        df=df,
        out_file=heatmap_file,
        top_terms=args.top_terms_heatmap,
        dpi=args.dpi,
    )

    print(f"KEGG figures written to: {out_dir}")


if __name__ == "__main__":
    main()
