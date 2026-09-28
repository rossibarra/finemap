# Changelog

## 0.4 — 2026-09-28

Fix the AGPv2 → v5 lift-over. The chain previously named `data/v2v5.chain`
actually mapped v5 → AGPv2, so every AGPv2 source was lifted in the wrong
direction. All maps and downstream analyses are rebuilt; the v4 → v5 lift-over
(Samayoa landrace and teosinte) was checked and is unaffected.

- Rename the AnchorWave chain to `data/v5v2.chain` and add
  `scripts/swap_chain.py`, which inverts it into a true AGPv2 → v5
  `data/v2v5.chain` (verified by an exact v5 → v2 → v5 round trip).
- Re-lift Rodgers-Melnick (88,863 → 118,339 intervals) and European HMM
  (21,026 → 27,328) crossovers; `data/jri_v5.bed` now has 409,525 intervals.
- Rebuild `data/finemap_v5.bed` (193,790 segments), `data/hapmap/`,
  `data/finemap_hierarchical_v5.bed`, and the hotspot BEDs (918 at 30×, 757 at
  1 kb/20×).
- Regenerate `data/ogut_v5.csv` (5,837 → 6,136 markers); it now agrees with an
  independent v2 → v4 → v5 lift of the Ogut map to within 10 kb for 98% of
  shared markers. `scripts/plot_marey_comparison.py` writes `Chr`-prefixed
  input for the new chain.
- Rerun the recombination, hotspot, π and gene-coverage analyses and update
  their numbers in the README and PI_*.md writeups. Correlations of π with
  recombination strengthen slightly; conclusions are unchanged apart from the
  wording noted in those documents.

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
