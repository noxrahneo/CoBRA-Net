#!/usr/bin/env python3

"""01_create_bigboss.py: building the BigBoss inventory from metadata (SampleStats.txt + Table EV4) and raw file names."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the BigBoss inventory from metadata (SampleStats.txt + Table EV4) and raw file names."
    )
    parser.add_argument(
        "--data-dir", default="data/GSE161529_RAW", help="Path to raw data folder with matrix/barcodes files."
    )
    parser.add_argument(
        "--sample-stats", default="data/HumanBreast10X-main/Tables/SampleStats.txt", help="Path to SampleStats.txt from the paper repo."
    )
    parser.add_argument(
        "--table-ev4", default="data/embj2020107333-sup-0006-tableev4.csv", help="Path to the Table EV4."
    )
    parser.add_argument(
        "--output", default="results/bigboss.csv", help="Output CSV path"
    )
    return parser.parse_args()

def load_and_validate_metadata(sample_stats_path: Path, table_ev4_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    #load both sources
    sample_stats = pd.read_csv(sample_stats_path, sep="\t")
    table_ev4 = pd.read_csv(table_ev4_path, skiprows=2)  #skip first two rows
    
    #validate sample names mathc in order
    print("Validating sample names")
    stats_names = sample_stats["SampleName"].tolist()
    ev4_names = table_ev4["Sample Name"].tolist()
    
    if len(stats_names) != len(ev4_names):
        raise ValueError(f"Sample name count mismatch: {len(stats_names)} in SampleStats vs {len(ev4_names)} in Table EV4")
    
    #check for each mathc
    mismatches = []
    for i, (s, e) in enumerate(zip(stats_names, ev4_names)):
        if s != e:
            mismatches.append((i, s, e))
    
    if mismatches:
        print(f"Found {len(mismatches)} sample name mismatches:")
        for i, s, e in mismatches:
            print(f"  Row {i}: SampleStats='{s}' vs TableEV4='{e}'")
        if len(mismatches) > 5:
            print(f"  ... and {len(mismatches) - 5} more")
    
    else:       print("All sample names match between SampleStats and Table EV4.")
    
    return sample_stats, table_ev4
    
def merge_metadata(sample_stats: pd.DataFrame, table_ev4: pd.DataFrame) -> pd.DataFrame:
    #merge on SampleName (avoid dublicating cell/gene counts)
    #drop dublicate columns from table_ev4
    
    ev4_clean = table_ev4.drop(columns=["Number of Cells", "Number of Cells After Filtrating", "Number of Genes Detected"], errors="ignore",)
    merged = pd.merge(sample_stats, ev4_clean, left_on="SampleName", right_on="Sample Name", how="left")
    
    return merged
     
def add_file_paths(bigboss: pd.DataFrame, data_dir: Path) -> pd.DataFrame:
    #find matrix/barcode files by GEO_ID
    rows_with_files = []
    for _, row in bigboss.iterrows():
        geo_id = row["GEO_ID"]
        
        mtx_files = list(data_dir.glob(f"{geo_id}*matrix.mtx.gz"))
        bc_files = list(data_dir.glob(f"{geo_id}*barcodes.tsv.gz"))
        
        #add a row for each matching file pair
        for mtx, bc in zip(mtx_files, bc_files):
            new_row = row.copy()
            new_row["MatrixFile"] = mtx.name
            new_row["BarcodesFile"] = bc.name
            rows_with_files.append(new_row)
    
    return pd.DataFrame(rows_with_files)
    
def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = (repo_root / args.data_dir).resolve()
    sample_stats_path = (repo_root / args.sample_stats).resolve()
    table_ev4_path = (repo_root / args.table_ev4).resolve()
    out_path = (repo_root / args.output).resolve()
    
    #validate file exists
    for path, name in [
        (sample_stats_path, "Sample Stats.txt"),
        (table_ev4_path, "Table EV4"),
        (data_dir, "Data Directory"),
    ]:
        if not path.exists():
            print(f"Error: {name} not found at {path}")
            return 1
    
    print (f"Loading metadata from both sources")
    sample_stats, table_ev4 = load_and_validate_metadata(sample_stats_path, table_ev4_path)
    
    print (f"Merging metadata")
    bigboss = merge_metadata(sample_stats, table_ev4)
    
    print("Adding file paths from raw data")
    bigboss = add_file_paths(bigboss, data_dir)
    
    #save
    out_path.parent.mkdir(parents=True, exist_ok=True)
    bigboss.to_csv(out_path, index=False)
    
    print(f"Saved: {len(bigboss)} entries to {out_path}")
    print(f"Rows: {len(bigboss)}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())