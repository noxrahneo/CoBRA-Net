#!/usr/bin/env python3
"""Tiny validation for Pearson correlation implementations.

Creates a toy matrix (10 samples x 5 genes) and compares:
1) Pairwise scipy.stats.pearsonr for all gene pairs
2) numpy.corrcoef(..., rowvar=False)
3) scipy.spatial.distance.pdist(..., metric="correlation") with 1 - distance

If all methods agree within tolerance, prints PASS.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr


def corr_by_pearsonr(x: np.ndarray) -> np.ndarray:
    """Compute full correlation matrix by explicit pearsonr pair loops."""
    n_genes = x.shape[1]
    corr = np.zeros((n_genes, n_genes), dtype=np.float64)
    for i in range(n_genes):
        corr[i, i] = 1.0
        for j in range(i + 1, n_genes):
            r, _ = pearsonr(x[:, i], x[:, j])
            corr[i, j] = r
            corr[j, i] = r
    return corr


def corr_by_numpy(x: np.ndarray) -> np.ndarray:
    corr = np.corrcoef(x, rowvar=False)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(corr, 1.0)
    return corr


def corr_by_scipy_distance(x: np.ndarray) -> np.ndarray:
    dist = pdist(x.T, metric="correlation")
    corr = 1.0 - squareform(dist)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(corr, 1.0)
    return corr


def main() -> None:
    rng = np.random.default_rng(42)

    # 10 samples x 5 genes
    x = rng.normal(size=(10, 5)).astype(np.float64)

    corr_pearsonr = corr_by_pearsonr(x)
    corr_numpy = corr_by_numpy(x)
    corr_scipy = corr_by_scipy_distance(x)

    diff_np = np.abs(corr_pearsonr - corr_numpy)
    diff_sp = np.abs(corr_pearsonr - corr_scipy)

    max_diff_np = float(diff_np.max())
    max_diff_sp = float(diff_sp.max())

    print("Toy shape (samples, genes):", x.shape)
    print("Max |pearsonr - np.corrcoef|:", max_diff_np)
    print("Max |pearsonr - scipy_distance|:", max_diff_sp)

    ok_np = np.allclose(corr_pearsonr, corr_numpy, atol=1e-12, rtol=1e-12)
    ok_sp = np.allclose(corr_pearsonr, corr_scipy, atol=1e-12, rtol=1e-12)

    if ok_np and ok_sp:
        print("PASS: All methods match within tolerance.")
    else:
        print("FAIL: Methods differ beyond tolerance.")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
