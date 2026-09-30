"""Quantify deformation between corresponding RMSD distance matrices."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


def validate_distance_matrix(matrix: np.ndarray, name: str) -> None:
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{name} must be a square matrix; got {matrix.shape}.")
    if not np.isfinite(matrix).all():
        raise ValueError(f"{name} contains NaN or infinite values.")
    if np.min(matrix) < -1e-10:
        raise ValueError(f"{name} contains negative distances.")
    if not np.allclose(matrix, matrix.T, atol=1e-7, rtol=0):
        raise ValueError(f"{name} is not symmetric within tolerance.")
    if not np.allclose(np.diag(matrix), 0.0, atol=1e-7, rtol=0):
        raise ValueError(f"{name} does not have a zero diagonal.")


def compare(gly: np.ndarray, ala: np.ndarray, resolution: int) -> pd.DataFrame:
    validate_distance_matrix(gly, "glycine matrix")
    validate_distance_matrix(ala, "alanine matrix")
    if gly.shape != ala.shape:
        raise ValueError(
            "The two matrices must have identical shape and matching torsional indexing. "
            f"Got {gly.shape} and {ala.shape}."
        )

    iu = np.triu_indices(gly.shape[0], 1)
    g = gly[iu].astype(float)
    a = ala[iu].astype(float)
    delta = a - g

    frob_g = np.linalg.norm(gly, "fro")
    delta_f = np.linalg.norm(ala - gly, "fro") / frob_g
    alpha = np.sum(gly * ala) / np.sum(gly * gly)
    delta_shape = np.linalg.norm(ala - alpha * gly, "fro") / frob_g

    values = {
        "resolution": f"{resolution}x{resolution}",
        "configurations": gly.shape[0],
        "pairwise_distances": g.size,
        "relative_frobenius_deformation": delta_f,
        "relative_frobenius_deformation_percent": 100 * delta_f,
        "optimal_global_scale": alpha,
        "global_scale_change_percent": 100 * (alpha - 1),
        "residual_shape_deformation": delta_shape,
        "residual_shape_deformation_percent": 100 * delta_shape,
        "pearson_correlation": pearsonr(g, a).statistic,
        "spearman_correlation": spearmanr(g, a).statistic,
        "contracted_pairs_percent": 100 * np.mean(delta < 0),
        "expanded_pairs_percent": 100 * np.mean(delta > 0),
        "unchanged_pairs_percent": 100 * np.mean(delta == 0),
        "mean_signed_change_angstrom": np.mean(delta),
        "mean_absolute_change_angstrom": np.mean(np.abs(delta)),
        "minimum_change_angstrom": np.min(delta),
        "maximum_change_angstrom": np.max(delta),
    }
    return pd.DataFrame([values])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gly", required=True, type=Path, help="Glycine RMSD .npy matrix")
    parser.add_argument("--ala", required=True, type=Path, help="Alanine RMSD .npy matrix")
    parser.add_argument("--resolution", required=True, type=int, choices=(25, 36))
    parser.add_argument("--output", required=True, type=Path, help="Output CSV path")
    args = parser.parse_args()

    gly = np.load(args.gly)
    ala = np.load(args.ala)
    result = compare(gly, ala, args.resolution)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(result.to_string(index=False))
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()

