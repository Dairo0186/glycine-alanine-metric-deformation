# Data availability

This repository contains the analysis code, common torsional index arrays, persistence intervals used for the 25 x 25 comparison, final figures and numerical summary tables.

The following large inputs are intentionally excluded:

- full 72 x 72 all-atom conformational ensembles;
- full 5184 x 5184 RMSD matrices;
- reduced 25 x 25 and 36 x 36 RMSD matrices.

They can be regenerated from the molecular-construction and RMSD workflow described in the manuscript. Every glycine and alanine matrix must use the same ordered torsional pairs. The ordering convention is

```text
q = 72 i + j
phi = 5 i degrees
psi = 5 j degrees
```

with `i,j = 0,...,71`. Thus, the first conformer is `(0 degrees, 0 degrees)` and the last is `(355 degrees, 355 degrees)`.

For a resolution `n x n`, common indices are selected using

```text
s_k = floor(72 k / n),  k = 0,...,n-1.
```

This correspondence is essential: equal matrix dimensions alone do not guarantee that rows and columns represent the same torsional conformers.

