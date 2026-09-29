# Review: outstanding concerns after the v2 → v5 lift-over fix (2026-09-28)

The old `data/v2v5.chain` actually converted v5 → AGPv2, so the AGPv2 crossover data (Rodgers-Melnick and European) and the Ogut markers were lifted in the wrong direction. The chain has been inverted (`scripts/swap_chain.py`, old file kept as `data/v5v2.chain`). `jri_v5.bed`, both maps, the hapmaps, the hotspots, `ogut_v5.csv`, and the downstream analyses were then rebuilt (see CHANGELOG 0.4). The items below are still open.

## Stale or unregenerated outputs

- **`results/finemap-map-comparison*.png` are stale.** They were made from the old maps, and no script in the repo generates them (they came from an HTML artifact). Regenerate or delete them.
- **`data/example1-4.bed` were not regenerated.** `scripts/simulate_example_regions.py` samples crossover interval lengths from `jri_v5.bed` at random, so these files still reflect the old interval set. The effect is probably minor, but they are out of date.
- **`scripts/plot_rate_along_chromosomes.py` was not rerun.** Its output (`results/rate_along_chromosomes.png`) is neither tracked nor referenced. It reads `ogut_v5.csv` and the maps, so rerun it if the figure is needed.

## Documentation / reproducibility

- **Resolved: the README's Step 2/3 commands didn't reproduce `jri_v5.bed`.** They wrote bare chromosome numbers for the Samayoa data, which `build_finemap.py` then silently skipped, and never assigned the `LRv4_`/`TEOv4_` IDs. `scripts/build_jri_v5.py` now does the whole conversion with validation and reproduces the tracked files. The map builders now stop with an error on a chromosome with no Ogut target.
- **Resolved: the recombination metaplot measured the wrong quantity.** `metaplot.py --uniform` divided each rate by its segment length, so values depended on the segmentation. It now uses an overlap-weighted mean, and the plot is regenerated (about 1.1–2.0 cM/Mb, previously 0.04–0.82).
- **Resolved: the `finemap` conda env was missing `openpyxl`**, even though `environment.yml` listed it. openpyxl 3.1.5 has now been installed into the env.
- **Three numbers in the PI writeups could not be reproduced, even from the old files:**
  - PI_VS_HOTSPOT_DISTANCE "3.1% of cM from intervals <100 bp"
  - PI_VS_HOTSPOT_DISTANCE masked/unmasked ρ "−0.755 / 1.000"
  - PI_CODING_VS_RECOMBINATION 0D "0.325" at ≥500 sites

  They were replaced with values from a re-implemented method. Confirm the method matches what was intended.

## Ogut map on v5

- **`data/ogut_v5.csv` and the independent `ogutweird/` hapmaps differ slightly.** `ogut_v5.csv` comes from a direct AnchorWave v2 → v5 lift; `ogutweird/` was built via v2 → v4 → v5.
  - **Agreement:** 5,683 markers are shared. Of these, 59% have identical positions, 98% are within 10 kb and 99.5% within 100 kb.
  - **Unshared markers:** 745 markers appear only in the hapmap and 453 only in `ogut_v5.csv`.
  - **Resolved: `ogut_v5.csv` is the correct one.** The check used the AGPv2 reference (`data/B73_RefGen_v2.fa.gz`, MaizeGDB, not tracked). For each marker, the 101 bp around its v2 position was searched for within ±300 bp of each candidate v5 position, on both strands.
    - These numbers predate the fix to `lift_ogut()` coordinates (see below) and could not be recomputed here because `ogutweird/` is not in the repository. The fix moves 36 reverse-strand markers by 2 bp and changes which markers lift (5 gained, 2 lost), so at most 36 of the shared markers can change class; the conclusion stands.
    - **Across all 5,683 shared markers:** the v2 sequence sits exactly at the `ogut_v5.csv` position for 97.3%, and at the ogutweird position for 58.7%.
    - **Where the two disagree (2,305 markers):** `ogut_v5.csv` alone is exact for 2,196, ogutweird alone for 5, both for 11 and neither for 93. The median mismatch at ogutweird's positions is 59 of 101 bp, i.e. no match.
    - **Conclusion:** the v2 → v4 → v5 route used for ogutweird misplaces about 40% of markers by bp to kb.
- **Resolved: 115 of 6,139 lifted Ogut markers (1.9%) lack sequence support at their position.** An earlier version reported 242 of 6,136 and called the 34 rejected markers that were out of cM order lift-over errors. That was wrong: `lift_ogut()` fed 1-based positions to CrossMap as BED starts and read the lifted BED start back as 1-based. The errors cancel on forward chains but put all 36 markers on reverse-strand chains 2 bp off (e.g. M395 at 66,452,031 instead of 66,452,033). With the fix, all 36 match exactly. Unsupported markers are dropped with `scripts/verify_ogut_v5.py --filter`, leaving 6,024.
  - **The ungapped ±50 bp test was too strict.** Of 209 markers that failed it, 90% sat within 50 bp of an alignment gap in the chain. The check now also accepts a marker when either one-sided flank (marker base plus 50 bp) matches exactly at the position with ≤2 mismatches. That recovers 94 markers (78 with 0 mismatches) and adds no cM-order breaks. The 115 still rejected have ≥3 mismatches on their better flank (median 23), a clean gap from the passing ones.
  - **35 markers are still out of cM order** (34 on reverse-strand chains plus M2180, next to M2178 on chr3). They lie in blocks that v5 inverts relative to AGPv2 (e.g. the last 2 Mb of chr2, 141–143 Mb on chr7 and 121.7–122.2 Mb on chr6). The Ogut cM values are monotone in AGPv2 position, so the Ogut order follows AGPv2 across these blocks. The markers are placed correctly by sequence and are kept. The inversions may reflect AGPv2 orientation errors; this has not been checked.
- **The corrected map fits the Ogut Marey curve worse on chr1, chr2 and chr7** (RMSE 2.6→2.8, 3.4→5.1 and 3.0→5.5 cM), though it fits better on the other seven. This is worth a look for local problems, e.g. rearrangements between AGPv2 and v5.

## Lift-over quality

- **The swapped chain is not netted from the v2 side.** An inverted v5 → v2 chain can contain overlapping blocks in v2 coordinates, so some v2 positions may map to more than one place. CrossMap's handling of those cases was not audited.
- **Some lifted intervals changed length a lot.** For 1% of re-lifted v2 intervals, v5 length / v2 length is above 2.96, and for 1% it is below 0.61. These are likely intervals spanning structural differences. Filtering them could sharpen the map.
- **About 13–16% of v2 intervals still fail to lift** (RM 87.0% retained, European 84.2%).
- **Resolved: crossover endpoints were split as `[start, start+1)` / `[end-1, end)` from 1-based source coordinates**, so every interval was 1 bp short. The left marker is now `[start-1, start)`. Result: 409,510 intervals (15 fewer, because endpoints next to alignment gaps flip between lifting and failing). The maps are identical at 1 Mb scale, but `finemap_v5.bed` now has 262,448 segments instead of 193,790, because an interval ending at a SNP and one starting at the same SNP now both include it, so their edges are 1 bp apart instead of shared, adding many 1-bp segments.

## Files not committed

- **`ogutweird/`** (external comparison hapmaps) and **`results/rate_chr4_83.5-87.5Mb.png`** were never tracked and are no longer on disk. The ogutweird comparison numbers above therefore can't be recomputed after the coordinate fix until those hapmaps are restored.
