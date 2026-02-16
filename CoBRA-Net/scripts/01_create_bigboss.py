#!/usr/bin/env python3

"""01_create_bigboss.py: building the BigBoss inventory from metadata and raw file names."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the BigBoss inventory from metadata and raw file names."
    )
    parser.add_argument(
        "--data-dir", default="data/GSE161529_RAW", help="Path to raw data folder with matrix/barcodes files."
    )
    parser.add_argument(
        "--metadata", default="data/HumanBreast10X-main/Tables/SampleStats.txt", help="Path to SamplesStats.txt from the paper repo."
    )
    parser.add_argument(
        "--output", default="results/bigboss.csv", help="Output CSV path"
    )
    return parser.parse_args()

def print_unique_titles (meta_path: Path) -> None:
    #quick check at the raw titles
    sample_stats = pd.read_csv(meta_path, sep="\t")
    titles = sorted(
        sample_stats["Title"].dropna().unique().tolist()
    )
    print(f"Unique title values({len(titles)}):")
    for t in titles:
        print(f"  - {t}")
    

def classify_sample_type(title: str) -> str: 
    title_lower = title.lower()
    
    if "normal" in title_lower:
        return "normal"
    if "brca1 pre-neoplastic" in title_lower:
        return "BRACA1_PreNeoplastic"
    if "triple negative brca1" in title_lower:
        return "TripleNegative_BRCA1"
    if "triple negative" in title_lower:
        return "TripleNegative"
    if "her2+" in title_lower:
        return "HER2_Positive"
    if "er+" in title_lower or "pr+" in title_lower:
        return "ER_Positive"
    
def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = (repo_root / args.data_dir).resolve()
    meta_path = (repo_root / args.metadata).resolve()
    out_path = (repo_root / args.output).resolve()
    
    #load metadata from the paper
    meta = pd.read_csv(meta_path, sep="\t")
    
    #match files by GEO_ID
    rows = []
    for _, row in meta.iterrows():
        geo_id = row["GEO_ID"]
        
        #find matching files in the data dir
        mtx_files = list(data_dir.glob(f"{geo_id}*matrix.mtx.gz"))
        bc_files = list(data_dir.glob(f"{geo_id}*barcodes.tsv.gz"))
        
        for mtx, bc in zip(mtx_files, bc_files):
            rows.append(
                {
                    "GEO_ID": geo_id,
                    "MatrixFile": mtx.name,
                    "BarcodesFile": bc.name,
                    "SampleName": row["SampleName"],
                    "Title": row["Title"],
                    "CellNumAfter": row["CellNumAfter"],
                    "GenesDetected": row["GenesDetected"],
                    "CellNum": row["CellNum"],
                    "Mito": row["Mito"],
                    "GeneLower": row["GeneLower"],
                    "GeneUpper": row["GeneUpper"],
                    "LibSize": row["LibSize"],
                }
            )
    bigboss = pd.DataFrame(rows)
    
    #add sample type
    bigboss["SampleType"] = bigboss["Title"].apply(classify_sample_type)
    
    #save
    out_path.parent.mkdir(parents=True, exist_ok=True)
    bigboss.to_csv(out_path, index=False)
    
    print(f"Saved: {len(bigboss)} entries to {out_path}")
    print(f"Rows: {len(bigboss)}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())