# Changelog

## 0.3 — 2026-09-01

New alternative recombination map: a smoothed, Ogut-calibrated hierarchical
interval-censored map, alongside the existing default interval-density map.

- Add `scripts/build_hierarchical_finemap.py`, which models each crossover as
  interval-censored under a shared chromosome-wide rate function, represents
  `log(r)` as piecewise constant in 100 kb bins, and applies a Gaussian
  random-walk penalty across neighboring bins. Fitted per chromosome by MAP
  optimization, then scaled to Ogut chromosome lengths. Options: `--bin-size`,
  `--smoothness`, `--iterations`, `--learning-rate`, `--output`.
- Add `data/finemap_hierarchical_v5.bed` — 21,325 non-overlapping 100 kb
  segments (shorter terminal segment per chromosome), columns `chrom`, `start`,
  `end`, `cM_start`, `cM_end`, `cM_per_Mb`.
- Add map-comparison plots in `results/`: genome-wide, chr8, chr10, and chr10
  140–145 Mb / 140–160 Mb zooms.
- `scripts/build_finemap.py`: raise `FloatingPointError` on decreasing cM
  instead of silently patching the interval and continuing.
- README: document both maps as Ogut-calibrated (crossover data set the
  within-chromosome profile, Ogut sets total cM per chromosome), add the
  "Alternative Hierarchical Map" section, with the interval-censored
  likelihood rendered in LaTeX.

Also included since 0.2 (June–July 2026):

- Add Ogut v5 genetic map data and `data/ogut_v5.csv`, exported from the plot
  script and documented in the README.
- Add a per-chromosome recombination rate plot as a standalone script.
- Replace the two-flag ogut-v5 regeneration path with a single
  `--regenerate-ogut-v5` flag; improve pipeline reproducibility docs and CLIs.

## 0.2 — 2026-04-17

- Add gene±1kb physical vs cM coverage analysis; `--flank` CLI arg (default
  1000 bp) on `scripts/plot_gene_cM_coverage.py`.
- README: gene±1kb section with plot and ±100 bp vs ±1 kb comparison table,
  plus a table of contents.

## 0.1

- Initial tagged release of the composite map pipeline.
