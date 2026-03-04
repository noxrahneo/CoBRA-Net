#!/usr/bin/env python3
"""Minimal cluster annotation for integrated condition objects."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc

from annotation_signature_utils import build_signatures
from h5ad_compat import read_h5ad_compat


REPO_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Minimal annotation: score signatures "
            "and map leiden clusters"
        )
    )
    parser.add_argument(
        "--condition",
        default="all",
        help="Condition name or 'all'",
    )
    parser.add_argument(
        "--list-conditions",
        action="store_true",
        help="List available condition folders and exit",
    )
    parser.add_argument(
        "--input-dir",
        default="results/stages/03_integration/integrated",
        help="Root with integrated condition folders",
    )
    parser.add_argument(
        "--output-dir",
        default="results/stages/04_annotation_minimal",
        help="Root directory for annotation outputs",
    )
    parser.add_argument(
        "--cluster-col",
        default="leiden",
        help="Cluster column in adata.obs",
    )
    parser.add_argument(
        "--signature-file",
        default="data/HumanBreast10X-main/Signatures/ImmuneMarkers2.txt",
        help="Optional TSV with columns CellType and Signatures",
    )
    parser.add_argument(
        "--lineage-rdata",
        default="data/HumanBreast10X-main/Signatures/Human-PosSigGenes.RData",
        help="RData file with lineage signatures",
    )
    parser.add_argument(
        "--disable-lineage-rdata",
        action="store_true",
        help="Disable lineage signatures from RData",
    )
    parser.add_argument(
        "--score-layer",
        choices=["X", "counts"],
        default="X",
        help="Matrix layer used by score_genes",
    )
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


def list_conditions(input_root: Path) -> list[str]:
    if not input_root.exists():
        return []
    out: list[str] = []
    for child in sorted(p for p in input_root.iterdir() if p.is_dir()):
        try:
            _ = integrated_file_for(child)
        except FileNotFoundError:
            continue
        out.append(child.name)
    return out


def resolve_conditions(input_root: Path, requested: str) -> list[str]:
    available = list_conditions(input_root)
    if not available:
        raise FileNotFoundError(
            "No integrated condition folders found in: "
            f"{input_root}"
        )
    if requested.strip().lower() == "all":
        return available
    if requested in available:
        return [requested]
    raise FileNotFoundError(
        f"Condition '{requested}' not found. Available: {', '.join(available)}"
    )


def score_signatures(
    adata: sc.AnnData,
    signatures: dict[str, list[str]],
    score_layer: str,
) -> list[str]:
    gene_pool = set(map(str, adata.var_names))
    score_cols: list[str] = []
    for label, genes in signatures.items():
        present = [gene for gene in genes if gene in gene_pool]
        if not present:
            continue
        score_col = f"sig_{label}"
        if score_layer == "counts" and "counts" in adata.layers:
            sc.tl.score_genes(
                adata,
                gene_list=present,
                score_name=score_col,
                use_raw=False,
                layer="counts",
            )
        else:
            sc.tl.score_genes(
                adata,
                gene_list=present,
                score_name=score_col,
                use_raw=False,
            )
        score_cols.append(score_col)
    return score_cols


def infer_mapping(
    adata: sc.AnnData,
    cluster_col: str,
    score_cols: list[str],
) -> pd.DataFrame:
    grouped = (
        adata.obs[[cluster_col, *score_cols]]
        .groupby(cluster_col, observed=False)
        .mean()
        .reset_index()
    )
    rows = []
    for _, row in grouped.iterrows():
        ordered = sorted(
            ((col, float(row[col])) for col in score_cols),
            key=lambda item: item[1],
            reverse=True,
        )
        top1_col, top1 = ordered[0]
        top2 = ordered[1][1] if len(ordered) > 1 else float("nan")
        rows.append(
            {
                cluster_col: str(row[cluster_col]),
                "predicted_cell_type": top1_col.replace("sig_", ""),
                "score_top1": top1,
                "score_top2": top2,
                "score_margin": (
                    top1 - top2 if pd.notna(top2) else float("nan")
                ),
            }
        )
    return pd.DataFrame(rows)


def annotate_one_condition(
    condition: str,
    input_root: Path,
    output_root: Path,
    signatures: dict[str, list[str]],
    cluster_col: str,
    score_layer: str,
) -> dict[str, object]:
    h5ad = integrated_file_for(input_root / condition)
    out_dir = output_root / condition
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = read_h5ad_compat(h5ad)
    if cluster_col not in adata.obs.columns:
        raise ValueError(
            f"Cluster column '{cluster_col}' missing in {condition}"
        )

    score_cols = score_signatures(adata, signatures, score_layer)
    if not score_cols:
        raise ValueError(
            f"No signature genes found in var_names for {condition}"
        )

    mapping_df = infer_mapping(adata, cluster_col, score_cols)
    mapping_df.to_csv(
        out_dir / f"{condition}_cluster_to_celltype_mapping.csv",
        index=False,
    )

    mapping = dict(
        zip(
            mapping_df[cluster_col].astype(str),
            mapping_df["predicted_cell_type"],
        )
    )
    adata.obs["cell_type_annot"] = (
        adata.obs[cluster_col].astype(str).map(mapping).fillna("Unknown")
    )

    adata.write_h5ad(out_dir / f"{condition}_annotated.h5ad")

    row = {
        "condition": condition,
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "clusters": int(adata.obs[cluster_col].nunique()),
        "annotated_cell_types": int(adata.obs["cell_type_annot"].nunique()),
        "integrated_input": str(h5ad),
        "annotated_output": str(out_dir / f"{condition}_annotated.h5ad"),
    }
    pd.DataFrame([row]).to_csv(
        out_dir / f"{condition}_annotation_summary.csv",
        index=False,
    )
    print(f"Annotation complete: {condition}")
    return row


def main() -> int:
    args = parse_args()
    input_root = resolve_base(args.input_dir)
    output_root = resolve_base(args.output_dir)
    output_root.mkdir(parents=True, exist_ok=True)

    if args.list_conditions:
        names = list_conditions(input_root)
        if not names:
            print("No conditions found")
            return 1
        print("Available conditions:")
        for name in names:
            print(f"- {name}")
        return 0

    conditions = resolve_conditions(input_root, args.condition)
    signatures = build_signatures(
        sig_path=resolve_base(args.signature_file),
        lineage_rdata_path=resolve_base(args.lineage_rdata),
        use_lineage_rdata=not args.disable_lineage_rdata,
        pam50_path=Path("/dev/null"),
        include_pam50=False,
        collapse_aliases=True,
        prune_redundant=True,
    )

    rows = []
    for condition in conditions:
        rows.append(
            annotate_one_condition(
                condition=condition,
                input_root=input_root,
                output_root=output_root,
                signatures=signatures,
                cluster_col=args.cluster_col,
                score_layer=args.score_layer,
            )
        )

    summary_file = output_root / "annotation_summary.csv"
    pd.DataFrame(rows).to_csv(summary_file, index=False)
    print(f"Summary written: {summary_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
