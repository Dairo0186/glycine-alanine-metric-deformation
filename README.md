# Glycine-Alanine Metric Deformation

Reproducible code and selected results for the manuscript:

**Metric Deformation and Topological Persistence in the Ideal Conformational Spaces of Ace-Gly-NMe and Ace-Ala-NMe**

## Scientific purpose

This repository quantifies how replacing a hydrogen at the alpha carbon of glycine with the methyl group of alanine changes the all-atom RMSD geometry of the corresponding ideal peptide configuration spaces.

Ace-Gly-NMe and Ace-Ala-NMe are sampled on the same periodic torsional domain defined by the backbone angles phi and psi. The comparison is performed between RMSD distance matrices indexed by the same torsional conformers, rather than by subtracting Cartesian vectors of unequal dimension.

The analysis separates:

- global metric deformation;
- uniform contraction or expansion;
- residual nonuniform deformation;
- local directional deformation and anisotropy;
- preservation or displacement of persistent-homology features.

## Main results

| Descriptor | 25 x 25 | 36 x 36 |
|---|---:|---:|
| Relative Frobenius deformation | 4.159% | 4.164% |
| Optimal global scale factor | 0.973616 | 0.973631 |
| Residual shape deformation | 3.215% | 3.223% |
| Pearson correlation | 0.993986 | 0.994001 |
| Spearman correlation | 0.990917 | 0.990892 |
| Contracted pairs in alanine | 74.002% | 73.909% |

For both molecules and both resolutions, the dominant persistent signature is

`(beta_0, beta_1, beta_2) = (1, 2, 1)`,

which is compatible with a two-dimensional torus. The methyl substitution therefore produces a small, nonuniform metric deformation while preserving the dominant topological organization within the studied resolutions.

## Repository structure

```text
.
├── data/
│   ├── indices/                 # Common deterministic torsional indices
│   └── persistence_25x25/       # Persistence intervals used in Figure 3
├── docs/                        # Notes on data availability
├── results/
│   ├── figures/                 # Final figures in PNG and PDF
│   └── tables/                  # Numerical summaries in CSV format
├── src/
│   ├── analyze_metric_deformation.py
│   └── generate_selected_figures.py
├── CITATION.cff
├── LICENSE
└── requirements.txt
```

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

## Metric analysis

The complete RMSD matrices are not stored in the repository because of their size. Place corresponding glycine and alanine matrices in a local directory and run:

```bash
python src/analyze_metric_deformation.py \
  --gly path/to/gly_RMSD_25x25.npy \
  --ala path/to/ala_RMSD_25x25.npy \
  --resolution 25 \
  --output results/tables/metric_summary_25x25_recomputed.csv
```

Repeat with the 36 x 36 matrices by changing the input paths and `--resolution 36`.

The script verifies shape equality, symmetry, finite entries, nonnegativity and zero diagonals before computing the comparison.

## Figures

To regenerate the persistence comparison from the selected intervals included here:

```bash
python src/generate_selected_figures.py
```

Figures 1 and 2 are distributed as final outputs. Their complete regeneration requires the RMSD matrices and full 72 x 72 local-deformation table described in `docs/DATA_AVAILABILITY.md`.

## Methodological antecedent

This study extends the correspondence-preserving framework introduced in:

Hernández, D. J., Cadavid, C. A., De Luque, J., Fernández Bueno, D., Vega, R. R., & Herrera, Á. R. (2026). Metric deformation and topological persistence of molecular configuration spaces. *Mathematical and Computational Applications, 31*(5), 185. https://doi.org/10.3390/mca31050185

The published work compares ideal and MMFF94-relaxed realizations of each molecular system. The present repository applies the same distinction between metric deformation and topological persistence to two different peptide models sharing a common torsional indexing.

## License

Code is released under the MIT License. Scientific results should be cited using the metadata in `CITATION.cff`.

