#!/usr/bin/env python3
"""Copy existing results into a staged layout (non-destructive)."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Copy results into results/stages (non-destructive)"
    )
    parser.add_argument(
        "--results-dir",
        default="results",
        help="Root results directory",
    )
    parser.add_argument(
        "--stages-dir",
        default="results/stages",
        help="Root staged results directory",
    )
    parser.add_argument(
        "--mode",
        choices=["copy", "symlink"],
        default="copy",
        help="Copy data or create symlinks",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print actions without writing",
    )
    return parser.parse_args()


def resolve(path_like: str) -> Path:
    path = Path(path_like)
    return path if path.is_absolute() else REPO_ROOT / path


def copy_tree(src: Path, dest: Path, dry_run: bool) -> None:
    if dry_run:
        print(f"[dry-run] copytree {src} -> {dest}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dest, dirs_exist_ok=True, copy_function=shutil.copy2)


def copy_file(src: Path, dest: Path, dry_run: bool) -> None:
    if dry_run:
        print(f"[dry-run] copy {src} -> {dest}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)


def link_path(src: Path, dest: Path, dry_run: bool) -> None:
    if dry_run:
        print(f"[dry-run] symlink {src} -> {dest}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        return
    dest.symlink_to(src, target_is_directory=src.is_dir())


def copy_annotation_conditions(
    src_root: Path,
    dest_root: Path,
    mode: str,
    dry_run: bool,
) -> None:
    for child in sorted(p for p in src_root.iterdir() if p.is_dir()):
        condition = child.name
        if (child / f"{condition}_annotated.h5ad").exists():
            if mode == "symlink":
                link_path(child, dest_root / condition, dry_run)
            else:
                copy_tree(child, dest_root / condition, dry_run)


def main() -> int:
    args = parse_args()
    results_dir = resolve(args.results_dir)
    stages_dir = resolve(args.stages_dir)

    mapping = [
        (
            results_dir / "r_to_python" / "qc_pre",
            stages_dir / "01_qc" / "qc_pre",
        ),
        (
            results_dir / "r_to_python" / "qc_post",
            stages_dir / "01_qc" / "qc_post",
        ),
        (
            results_dir / "filtered_samples",
            stages_dir / "02_preprocess" / "filtered",
        ),
        (
            results_dir / "preprocessed_samples",
            stages_dir / "02_preprocess" / "preprocessed",
        ),
        (
            results_dir / "integrated_samples",
            stages_dir / "03_integration" / "integrated",
        ),
        (
            results_dir / "annotation" / "composition",
            stages_dir / "05_composition",
        ),
        (
            results_dir / "annotation" / "kegg",
            stages_dir / "06_kegg",
        ),
        (
            results_dir / "annotation" / "thesis_pack",
            stages_dir / "90_reports" / "thesis_pack",
        ),
    ]

    for src, dest in mapping:
        if not src.exists():
            print(f"[skip] missing {src}")
            continue
        if args.mode == "symlink":
            link_path(src, dest, args.dry_run)
        else:
            copy_tree(src, dest, args.dry_run)

    # Annotation stage: copy condition folders only.
    annot_src = results_dir / "annotation"
    annot_dest = stages_dir / "04_annotation"
    if annot_src.exists():
        copy_annotation_conditions(
            annot_src, annot_dest, args.mode, args.dry_run
        )

    # Reports: copy high-level summary files if present.
    for name in ["marker_heatmap_summary.csv", "thesis_one_page_summary.md"]:
        src = results_dir / "annotation" / name
        if src.exists():
            dest = stages_dir / "90_reports" / name
            if args.mode == "symlink":
                link_path(src, dest, args.dry_run)
            else:
                copy_file(src, dest, args.dry_run)

    print("Migration complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
