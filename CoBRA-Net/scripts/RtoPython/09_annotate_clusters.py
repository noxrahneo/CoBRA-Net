#!/usr/bin/env python3
"""Annotate integrated clusters for one or all conditions.

R_TRANSLATED: yes

R -> Python mapping used here:
- R: FindMarkers per cluster
    Python: sc.tl.rank_genes_groups + top marker export.
- R: marker/signature heatmap checks (ImmuneMarkers2, curated sets)
    Python: sc.tl.score_genes + cluster-level signature means.
- R: cluster interpretation from markers/signatures
    Python: cluster -> cell-type mapping table + marker-overlap sanity check.
- R: tSNE by cluster/group figures
    Python: UMAP regenerated with predicted cell-type labels.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "analysis"))

from warehouse import (  # noqa: E402
    WarehouseRecord,
    append_warehouse,
    params_hash,
    utc_now_iso,
)

DEFAULT_SIGNATURES: dict[str, list[str]] = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "MSLN"],
    "Basal_Epi": ["KRT5", "KRT14", "TP63", "ITGA6", "KRT17"],
    "Luminal_Epi": ["KRT8", "KRT18", "ESR1", "PGR", "GATA3"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "CLDN5"],
    "TCell": ["CD3D", "CD3E", "TRAC", "IL7R", "LTB"],
    "NK": ["NKG7", "GNLY", "KLRD1", "TRBC1", "TRBC2"],
    "BCell": ["MS4A1", "CD79A", "CD74", "CD19", "HLA-DRA"],
    "Plasma": ["MZB1", "JCHAIN", "SDC1", "XBP1", "DERL3"],
    "Myeloid": ["LYZ", "LST1", "FCER1G", "CTSS", "TYROBP"],
    "Mast": ["TPSAB1", "TPSB2", "KIT", "MS4A2", "HDC"],
    "Cycling": ["MKI67", "TOP2A", "PCNA", "CENPF", "TYMS"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Annotate integrated clusters for one "
            "condition or all conditions"
        )
    )
    parser.add_argument(
        "--condition",
        default="Normal",
        help=(
            "Condition folder name under "
            "results/stages/03_integration/integrated, or 'all'"
        ),
    )
    parser.add_argument(
        "--list-conditions",
        action="store_true",
        help="List available condition folders and exit",
    )
    parser.add_argument(
        "--input-dir",
        default="results/stages/03_integration/integrated",
        help="Root directory with integrated condition outputs",
    )
    parser.add_argument(
        "--output-dir",
        default="results/stages/04_annotation",
        help="Directory to write annotation outputs",
    )
    parser.add_argument(
        "--cluster-col",
        default="leiden",
        help="Cluster column in adata.obs",
    )
    parser.add_argument(
        "--n-top-markers",
        type=int,
        default=25,
        help="Top DE markers per cluster to export",
    )
    parser.add_argument(
        "--signature-file",
        default="data/HumanBreast10X-main/Signatures/ImmuneMarkers2.txt",
        help="Optional TSV with columns CellType and Signatures",
    )
    parser.add_argument(
        "--score-layer",
        choices=["X", "counts"],
        default="X",
        help="Matrix source for score_genes",
    )
    parser.add_argument(
        "--previous-label-col",
        default="",
        help="Optional existing label column to compare with annotation",
    )
    parser.add_argument(
        "--title-prefix",
        default="",
        help="Optional custom plot title prefix",
    )
    parser.add_argument("--dpi", type=int, default=160, help="Figure DPI")
    return parser.parse_args()


def resolve_base(path_like: str) -> Path:
    path = Path(path_like)
    return path if path.is_absolute() else REPO_ROOT / path


def integrated_file_for(condition_dir: Path) -> Path:
    preferred = condition_dir / f"{condition_dir.name}_integrated.h5ad"
    if preferred.exists():
        return preferred
    hits = sorted(condition_dir.glob("*_integrated.h5ad"))
    if not hits:
        raise FileNotFoundError(
            f"No integrated h5ad found in: {condition_dir}"
        )
    return hits[0]


def list_conditions(in_root: Path) -> list[str]:
    if not in_root.exists():
        return []
    names: list[str] = []
    for child in sorted(p for p in in_root.iterdir() if p.is_dir()):
        try:
            _ = integrated_file_for(child)
        except FileNotFoundError:
            continue
        names.append(child.name)
    return names


def resolve_conditions(in_root: Path, requested: str) -> list[str]:
    available = list_conditions(in_root)
    if not available:
        raise FileNotFoundError(
            f"No condition folders with integrated files found in: {in_root}"
        )

    query = requested.strip()
    if query.lower() == "all":
        return available
    if query in available:
        return [query]
    raise FileNotFoundError(
        f"Condition '{requested}' not found. Available: {', '.join(available)}"
    )


def load_signature_table(path: Path) -> dict[str, list[str]]:
    """Read optional external signature table from publication resources.

    Expected columns: CellType, Signatures (one gene per row).
    """
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t")
    need = {"CellType", "Signatures"}
    if not need.issubset(set(df.columns)):
        return {}

    out: dict[str, list[str]] = {}
    for celltype, sub in df.groupby("CellType"):
        genes = (
            sub["Signatures"]
            .astype(str)
            .str.strip()
            .replace("", np.nan)
            .dropna()
            .tolist()
        )
        if genes:
            out[str(celltype)] = list(dict.fromkeys(genes))
    return out


def build_signatures(sig_path: Path) -> dict[str, list[str]]:
    """Merge built-in core signatures with optional external signatures.

    R logic equivalent: combine curated marker panels before cluster
    interpretation.
    """
    external = load_signature_table(sig_path)
    signatures = {k: list(v) for k, v in DEFAULT_SIGNATURES.items()}
    for key, genes in external.items():
        if key in signatures:
            signatures[key] = list(dict.fromkeys([*signatures[key], *genes]))
        else:
            signatures[key] = list(genes)
    return signatures


def ensure_annotation_inputs(adata: sc.AnnData, cluster_col: str) -> None:
    """Fail fast if required clustering/embedding fields are missing."""
    if cluster_col not in adata.obs.columns:
        raise ValueError(
            f"Cluster column '{cluster_col}' not found in obs"
        )
    if "X_umap" not in adata.obsm:
        raise ValueError(
            "X_umap missing. Run integration plotting step first."
        )


def maybe_write_confusion(
    adata: sc.AnnData,
    previous_label_col: str,
    out_file: Path,
) -> None:
    """Optional validation table against a previous cell label column."""
    if not previous_label_col or previous_label_col not in adata.obs.columns:
        return
    confusion = pd.crosstab(
        adata.obs[previous_label_col].astype(str),
        adata.obs["cell_type_annot"].astype(str),
        dropna=False,
    )
    confusion.to_csv(out_file)


def run_markers(
    adata: sc.AnnData,
    cluster_col: str,
    n_top_markers: int,
) -> pd.DataFrame:
    # R equivalent: FindMarkers per cluster (positive markers support).
    # Important: run DE on expression-like values,
    # not scaled integrated values.
    # If counts layer exists, rebuild log-normalized matrix for DE testing.
    adata_de = adata.copy()
    if "counts" in adata_de.layers:
        adata_de.X = adata_de.layers["counts"].copy()
        sc.pp.normalize_total(adata_de, target_sum=1e4)
        sc.pp.log1p(adata_de)

    sc.tl.rank_genes_groups(
        adata_de,
        groupby=cluster_col,
        method="wilcoxon",
        use_raw=False,
        n_genes=max(50, n_top_markers),
    )
    df = sc.get.rank_genes_groups_df(adata_de, group=None)
    df = df.rename(columns={"group": "cluster"})
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.sort_values(["cluster", "pvals_adj", "logfoldchanges"])
    top = df.groupby("cluster", as_index=False).head(n_top_markers)
    return top


def score_signatures(
    adata: sc.AnnData,
    signatures: dict[str, list[str]],
    score_layer: str,
) -> tuple[list[str], dict[str, int]]:
    # Signature-scoring support for cluster interpretation.
    # This mirrors the R workflow's marker/signature validation logic.
    if score_layer == "counts" and "counts" in adata.layers:
        source = adata.layers["counts"].copy()
        adata.layers["_score_tmp"] = source
        use_layer = "_score_tmp"
    else:
        use_layer = None

    score_cols: list[str] = []
    genes_used: dict[str, int] = {}
    gene_pool = set(map(str, adata.var_names))
    for label, genes in signatures.items():
        present = [g for g in genes if g in gene_pool]
        genes_used[label] = len(present)
        if not present:
            continue
        score_col = f"sig_{label}"
        sc.tl.score_genes(
            adata,
            gene_list=present,
            score_name=score_col,
            use_raw=False,
            layer=use_layer,
        )
        score_cols.append(score_col)

    if "_score_tmp" in adata.layers:
        del adata.layers["_score_tmp"]

    return score_cols, genes_used


def cluster_signature_table(
    adata: sc.AnnData,
    cluster_col: str,
    score_cols: list[str],
) -> pd.DataFrame:
    if not score_cols:
        return pd.DataFrame(columns=[cluster_col])
    table = (
        adata.obs[[cluster_col, *score_cols]]
        .groupby(cluster_col, observed=False)
        .mean()
        .reset_index()
    )
    return table


def infer_cluster_labels(
    score_table: pd.DataFrame,
    cluster_col: str,
) -> pd.DataFrame:
    score_cols = [c for c in score_table.columns if c.startswith("sig_")]
    if not score_cols:
        return pd.DataFrame(columns=[cluster_col, "predicted_cell_type"])

    records = []
    for _, row in score_table.iterrows():
        vals = row[score_cols].sort_values(ascending=False)
        top1 = vals.index[0]
        top2 = vals.index[1] if len(vals) > 1 else top1
        records.append(
            {
                cluster_col: str(row[cluster_col]),
                "predicted_cell_type": top1.replace("sig_", ""),
                "score_top1": float(vals.iloc[0]),
                "score_top2": float(vals.iloc[1]) if len(vals) > 1 else np.nan,
                "score_margin": float(vals.iloc[0] - vals.iloc[1])
                if len(vals) > 1
                else np.nan,
                "top1_signature": top1,
                "top2_signature": top2,
            }
        )
    return pd.DataFrame(records)


def marker_sanity(
    mapping_df: pd.DataFrame,
    marker_df: pd.DataFrame,
    signatures: dict[str, list[str]],
    cluster_col: str,
) -> pd.DataFrame:
    # Quick sanity check: overlap between inferred label signature
    # and top DE markers in the same cluster.
    top_markers = (
        marker_df.groupby("cluster", observed=False)["names"]
        .apply(list)
        .to_dict()
    )
    rows = []
    for _, row in mapping_df.iterrows():
        cl = str(row[cluster_col])
        label = str(row["predicted_cell_type"])
        sig = signatures.get(label, [])
        markers = top_markers.get(cl, [])
        overlap = sorted(set(markers).intersection(sig))
        rows.append(
            {
                cluster_col: cl,
                "predicted_cell_type": label,
                "marker_overlap_n": len(overlap),
                "marker_overlap_genes": ";".join(overlap),
            }
        )
    return pd.DataFrame(rows)


def make_umap_plot(
    adata: sc.AnnData,
    color: str,
    title: str,
    out_file: Path,
    dpi: int,
) -> None:
    fig = sc.pl.umap(
        adata,
        color=color,
        show=False,
        frameon=False,
        return_fig=True,
        title=title,
        legend_loc="right margin",
        legend_fontsize=8,
    )
    fig.savefig(
        out_file,
        dpi=dpi,
        bbox_inches="tight",
        pad_inches=0.1,
    )
    plt.close(fig)


def annotate_one_condition(
    args: argparse.Namespace,
    in_root: Path,
    out_root: Path,
    sig_path: Path,
    condition: str,
) -> dict[str, object]:
    """Run the full annotation workflow for a single integrated condition."""
    cond_dir = in_root / condition
    h5ad = integrated_file_for(cond_dir)

    out_dir = out_root / condition
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(h5ad)
    ensure_annotation_inputs(adata, args.cluster_col)

    # Step 1 (R FindMarkers equivalent): cluster marker discovery.
    marker_df = run_markers(
        adata=adata,
        cluster_col=args.cluster_col,
        n_top_markers=args.n_top_markers,
    )
    marker_df.to_csv(
        out_dir / f"{condition}_cluster_markers_top.csv",
        index=False,
    )

    # Step 2: signature table assembly (built-in + optional external).
    signatures = build_signatures(sig_path)

    # Step 3: score signatures and aggregate at cluster level.
    score_cols, genes_used = score_signatures(
        adata=adata,
        signatures=signatures,
        score_layer=args.score_layer,
    )

    score_table = cluster_signature_table(
        adata=adata,
        cluster_col=args.cluster_col,
        score_cols=score_cols,
    )
    score_table.to_csv(
        out_dir / f"{condition}_cluster_signature_scores.csv",
        index=False,
    )

    # Step 4: infer cluster -> cell-type mapping from signature maxima.
    mapping_df = infer_cluster_labels(
        score_table=score_table,
        cluster_col=args.cluster_col,
    )

    if mapping_df.empty:
        raise ValueError("No signature scores were computed for annotation.")

    # Step 5: marker-overlap sanity support, then persist mapping.
    sanity_df = marker_sanity(
        mapping_df=mapping_df,
        marker_df=marker_df,
        signatures=signatures,
        cluster_col=args.cluster_col,
    )

    mapping_df = mapping_df.merge(
        sanity_df,
        on=[args.cluster_col, "predicted_cell_type"],
        how="left",
    )
    mapping_df.to_csv(
        out_dir / f"{condition}_cluster_to_celltype_mapping.csv",
        index=False,
    )

    # Apply cluster label mapping at single-cell level.
    cluster_to_label = dict(
        zip(
            mapping_df[args.cluster_col].astype(str),
            mapping_df["predicted_cell_type"],
        )
    )
    adata.obs["cell_type_annot"] = (
        adata.obs[args.cluster_col]
        .astype(str)
        .map(cluster_to_label)
        .fillna("Unknown")
    )

    used_df = pd.DataFrame(
        {
            "signature": list(genes_used.keys()),
            "genes_present_n": list(genes_used.values()),
        }
    ).sort_values(["genes_present_n", "signature"], ascending=[False, True])
    used_df.to_csv(
        out_dir / f"{condition}_signature_gene_coverage.csv",
        index=False,
    )

    # Step 6: annotated UMAP exports (cell type + cluster + sample).
    title_prefix = args.title_prefix.strip() or condition
    make_umap_plot(
        adata,
        color="cell_type_annot",
        title=f"{title_prefix} | UMAP: cell_type_annot",
        out_file=fig_dir / f"{condition}_umap_cell_type_annot.png",
        dpi=args.dpi,
    )
    make_umap_plot(
        adata,
        color=args.cluster_col,
        title=f"{title_prefix} | UMAP: {args.cluster_col}",
        out_file=fig_dir / f"{condition}_umap_{args.cluster_col}.png",
        dpi=args.dpi,
    )

    if "SampleName" in adata.obs.columns:
        make_umap_plot(
            adata,
            color="SampleName",
            title=f"{title_prefix} | UMAP: SampleName",
            out_file=fig_dir / f"{condition}_umap_sample.png",
            dpi=args.dpi,
        )

    maybe_write_confusion(
        adata=adata,
        previous_label_col=args.previous_label_col,
        out_file=out_dir / f"{condition}_confusion_previous_vs_annot.csv",
    )

    # Step 7: save annotated object + compact summary row.
    log1p = adata.uns.get("log1p")
    if isinstance(log1p, dict) and "base" in log1p:
        try:
            log1p["base"] = float(log1p["base"])
        except (TypeError, ValueError):
            log1p["base"] = float(np.e)
    annotated_file = out_dir / f"{condition}_annotated.h5ad"
    adata.write_h5ad(annotated_file)

    summary_row = {
        "condition": condition,
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "clusters": int(adata.obs[args.cluster_col].nunique()),
        "annotated_cell_types": int(
            adata.obs["cell_type_annot"].nunique()
        ),
        "integrated_input": str(h5ad),
        "annotated_output": str(annotated_file),
    }
    pd.DataFrame([summary_row]).to_csv(
        out_dir / f"{condition}_annotation_summary.csv",
        index=False,
    )

    print(f"Annotation complete: {condition}")
    print(f"Output directory: {out_dir}")
    print(f"Annotated object: {annotated_file}")
    return summary_row


def main() -> int:
    args = parse_args()
    in_root = resolve_base(args.input_dir)
    out_root = resolve_base(args.output_dir)
    sig_path = resolve_base(args.signature_file)
    available = list_conditions(in_root)
    if args.list_conditions:
        if not available:
            print("No conditions found")
            return 0
        print("Available conditions:")
        for name in available:
            print(f"- {name}")
        return 0

    conditions = resolve_conditions(in_root, args.condition)
    arg_hash = params_hash(vars(args))
    now = utc_now_iso()
    rows: list[dict[str, object]] = []
    warehouse_rows: list[WarehouseRecord] = []
    for condition in conditions:
        row = annotate_one_condition(
            args=args,
            in_root=in_root,
            out_root=out_root,
            sig_path=sig_path,
            condition=condition,
        )
        rows.append(row)
        warehouse_rows.append(
            WarehouseRecord(
                input_file=str(row["integrated_input"]),
                output_file=str(row["annotated_output"]),
                script="scripts/RtoPython/09_annotate_clusters.py",
                date_utc=now,
                params_hash=arg_hash,
                condition=str(row["condition"]),
                stage="annotation",
            )
        )

    if len(rows) > 1:
        summary_all = pd.DataFrame(rows)
        summary_file = out_root / "annotation_summary.csv"
        summary_all.to_csv(summary_file, index=False)
        print(f"Combined summary: {summary_file}")

    warehouse_file = append_warehouse(out_root, warehouse_rows)
    print(f"Warehouse log: {warehouse_file}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
