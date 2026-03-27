# quinn_metagene.py changelog

This file summarizes how the current `quinn_metagene.py` differs from the original committed baseline.

## Major differences from the original

- the hard-coded gene BED path was removed
- gene annotation is now supplied at runtime with either `--gff3` or `--gene-bed`
- positional arguments were replaced with explicit command-line flags
- `--dsb-bed` is now optional, so the script can run as a single-track metagene plot
- signal BED parsing is more tolerant:
  - 3-column BED files assign a count of `1` per row
  - 4-column BED files with a non-numeric fourth column also assign a count of `1` for that row
- chromosome names are normalized before matching genes and signal tracks
- aggregation was refactored from per-gene pandas masking to chromosome-grouped NumPy arrays with sorted midpoint lookups for better speed
- output values are now normalized to density per bp so bins with different widths are comparable
- the output path is configurable with `--output`
- plotting can run without opening a window by using `--no-show`

## Algorithm and layout

The metagene layout itself is still the same:
- upstream flank
- TSS inner window
- scaled gene body
- TTS inner window
- downstream flank

The midpoint-based assignment of BED intervals to metagene bins is also unchanged. The main computational changes are performance-oriented refactoring and bin-width normalization of the plotted values.

## Original baseline behavior

The original script:
- required `co.bed dsb.bed 500 50 100` as positional arguments
- always read genes from `/home2/qyj2/from_maizegdb/v5/b73_v5_gene_tss_tts_sort.bed`
- required both CO and DSB inputs
- assumed a numeric count column in the BED inputs
- reported raw aggregated counts rather than density normalized by bin width
