"""Regenerate the selected persistence-diagram comparison included in the repository."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "persistence_25x25"
OUT = ROOT / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

COLORS = {0: "#0072B2", 1: "#D55E00", 2: "#009E73"}


def main() -> None:
    diagrams = {
        mol: [np.load(DATA / f"{mol}_PH_25x25_H{dim}.npy") for dim in range(3)]
        for mol in ("gly", "ala")
    }
    maximum = max(
        np.nanmax(array[np.isfinite(array)])
        for molecular_diagrams in diagrams.values()
        for array in molecular_diagrams
    ) * 1.04

    fig, axes = plt.subplots(1, 2, figsize=(9.3, 4.2), constrained_layout=True)
    for ax, mol, title in zip(axes, ("gly", "ala"), ("Glycine", "Alanine")):
        ax.plot([0, maximum], [0, maximum], "--", color="0.45", linewidth=1)
        for dim, intervals in enumerate(diagrams[mol]):
            finite = intervals[np.isfinite(intervals[:, 1])]
            ax.scatter(
                finite[:, 0], finite[:, 1], s=16 if dim else 10, alpha=0.72,
                color=COLORS[dim], edgecolors="none", label=rf"$H_{dim}$"
            )
        ax.set(
            xlim=(0, maximum), ylim=(0, maximum),
            xlabel="Birth (angstrom)", ylabel="Death (angstrom)",
            title=rf"{title}: 25$\times$25",
        )
        ax.set_aspect("equal", adjustable="box")
        ax.legend(loc="lower right", fontsize=8)

    for extension in ("png", "pdf"):
        path = OUT / f"figure_3_persistence_comparison_25x25.{extension}"
        fig.savefig(path, dpi=400 if extension == "png" else None, bbox_inches="tight")
        print(f"Saved: {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()

