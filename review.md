# Review: outstanding concerns after the v2 → v5 lift-over fix (2026-09-28)

The old `data/v2v5.chain` actually converted v5 → AGPv2, so the AGPv2 crossover data (Rodgers-Melnick and European) and the Ogut markers were lifted in the wrong direction. The chain has been inverted (`scripts/swap_chain.py`, old file kept as `data/v5v2.chain`). `jri_v5.bed`, both maps, the hapmaps, the hotspots, `ogut_v5.csv`, and the downstream analyses were then rebuilt (see CHANGELOG 0.4). The items below are still open.

## Stale or unregenerated outputs

- **`results/finemap-map-comparison*.png` are stale.** They were made from the old maps, and no script in the repo generates them (they came from an HTML artifact). Regenerate or delete them.
- **`data/example1-4.bed` were not regenerated.** `scripts/simulate_example_regions.py` samples crossover interval lengths from `jri_v5.bed` at random, so these files still reflect the old interval set. The effect is probably minor, but they are out of date.
- **`scripts/plot_rate_along_chromosomes.py` was not rerun.** Its output (`results/rate_along_chromosomes.png`) is neither tracked nor referenced. It reads `ogut_v5.csv` and the maps, so rerun it if the figure is needed.

## Documentation / reproducibility

- **The README's Step 3 command doesn't reproduce `jri_v5.bed`.** It concatenates the Samayoa `_v5.bed` files, which have 4 columns and no ID, but `jri_v5.bed` has 5 columns with `LRv4_`/`TEOv4_` IDs. The script or command that assigned those IDs is not documented. For this rebuild, the v4 rows were kept verbatim from the previous `jri_v5.bed`. A fresh CrossMap lift reproduced their coordinates exactly (263,858 of 263,858).
- **The `finemap` conda env lacks `openpyxl`**, even though `environment.yml` lists it. Step 1 (`scripts/hmm_co_pipeline.py`) was run with base Anaconda Python instead. Rebuild or update the env.
- **Three numbers in the PI writeups could not be reproduced, even from the old files:**
  - PI_VS_HOTSPOT_DISTANCE "3.1% of cM from intervals <100 bp"
  - PI_VS_HOTSPOT_DISTANCE masked/unmasked ρ "−0.755 / 1.000"
  - PI_CODING_VS_RECOMBINATION 0D "0.325" at ≥500 sites

  They were replaced with values from a re-implemented method. Confirm the method matches what was intended.

## Ogut map on v5

- **`data/ogut_v5.csv` and the independent `ogutweird/` hapmaps differ slightly.** `ogut_v5.csv` comes from a direct AnchorWave v2 → v5 lift; `ogutweird/` was built via v2 → v4 → v5.
  - **Agreement:** 5,683 markers are shared. Of these, 59% have identical positions, 98% are within 10 kb and 99.5% within 100 kb.
  - **Unshared markers:** 745 markers appear only in the hapmap and 453 only in `ogut_v5.csv`.
  - **Can't tell which is closer:** every disagreeing marker falls between the same flanking agreed markers under both liftovers. A base-level check would need the AGPv2 reference.
- **35 markers in `ogut_v5.csv` are out of cM order along v5.** The hapmap is fully ordered, so they were probably removed there. Consider filtering them.
- **The corrected map fits the Ogut Marey curve worse on chr1, chr2 and chr7** (RMSE 2.6→2.8, 3.4→5.1 and 3.0→5.5 cM), though it fits better on the other seven. This is worth a look for local problems, e.g. rearrangements between AGPv2 and v5.

## Lift-over quality

- **The swapped chain is not netted from the v2 side.** An inverted v5 → v2 chain can contain overlapping blocks in v2 coordinates, so some v2 positions may map to more than one place. CrossMap's handling of those cases was not audited.
- **Some lifted intervals changed length a lot.** For 1% of re-lifted v2 intervals, v5 length / v2 length is above 2.96, and for 1% it is below 0.61. These are likely intervals spanning structural differences. Filtering them could sharpen the map.
- **About 13–16% of v2 intervals still fail to lift** (RM 87.0% retained, European 84.2%).

## Files not committed

- **`ogutweird/`** — external comparison hapmaps.
- **`results/rate_chr4_83.5-87.5Mb.png`** — no script generates it.
