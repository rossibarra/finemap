# Changelog

## 0.4 — 2026-09-28

Fix the AGPv2 → v5 lift-over. The chain previously named `data/v2v5.chain`
actually mapped v5 → AGPv2, so every AGPv2 source was lifted in the wrong
direction. All maps and downstream analyses are rebuilt; the v4 → v5 lift-over
(Samayoa landrace and teosinte) was checked and is unaffected.

- Rename the AnchorWave chain to `data/v5v2.chain` and add
  `scripts/swap_chain.py`, which inverts it into a true AGPv2 → v5
  `data/v2v5.chain` (verified by an exact v5 → v2 → v5 round trip).
- Re-lift Rodgers-Melnick (88,863 → 118,323 intervals) and European HMM
  (21,026 → 27,328) crossovers; `data/jri_v5.bed` now has 409,510 intervals.
- Rebuild `data/finemap_v5.bed` (262,448 segments), `data/hapmap/`,
  `data/finemap_hierarchical_v5.bed`, and the hotspot BEDs (916 at 30×, 755 at
  1 kb/20×).
- Regenerate `data/ogut_v5.csv` and add `scripts/verify_ogut_v5.py`, which
  checks each lifted marker against the AGPv2 and v5 sequence. Markers without
  sequence support at their v5 position are dropped (6,139 lifted → 6,024 kept). A marker
  passes on the full ±50 bp window or on either one-sided flank (≤2 mismatches),
  so markers next to an assembly indel are not rejected.
  `lift_ogut()` now passes 1-based positions to CrossMap as `[position−1,
  position)` and adds 1 to the lifted start; before, reverse-strand chains put
  36 markers 2 bp off and the check wrongly rejected them. The 35 markers out
  of cM order are kept: they lie in blocks inverted between AGPv2 and v5. The sequence check also
  shows that an independent v2 → v4 → v5 lift of the Ogut map misplaces about
  40% of markers. `scripts/plot_marey_comparison.py` writes `Chr`-prefixed
  input for the new chain.
- Add `scripts/build_jri_v5.py`, which replaces the README's awk steps for
  the lift-over and the `jri_v5.bed` build. It treats source coordinates as 1-based
  (left endpoint marker `[start-1, start)`; the old awk used `[start, start+1)`,
  so every interval was 1 bp short) and validates chromosome names, schema, IDs and per-source counts.
  The old README commands wrote bare chromosome numbers for the Samayoa data and
  never assigned the `LRv4_`/`TEOv4_` IDs. `build_finemap.py` and
  `build_hierarchical_finemap.py` now stop with an error when a chromosome has
  no Ogut target, instead of silently skipping it.
- Fix `scripts/metaplot.py --uniform`, which divided each rate by its segment
  length, so bin values depended on how the map was segmented. Rates are now
  averaged over each bin weighted by overlap, dividing by the bp the map covers.
  The old behaviour, suited to per-interval totals, is kept as
  `--uniform-mode total`, and `--self-test` checks both. The recombination
  metaplot now reads about 1.1–2.0 cM/Mb (was 0.04–0.82) with a flatter
  gene-versus-flank contrast.
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
