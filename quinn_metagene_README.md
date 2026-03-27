# quinn_metagene.py

`quinn_metagene.py` plots CO and DSB dual-anchor metagene profiles.

It keeps the original profile layout:
- upstream flank
- TSS inner window
- scaled gene body
- TTS inner window
- downstream flank

## What changed

- the hard-coded gene BED path was removed
- gene annotations can now be read directly from `GFF3` with `--gff3`
- the original BED-style gene file is still supported with `--gene-bed`
- positional arguments were replaced with explicit command-line flags
- the aggregation code now uses chromosome-sorted arrays and interval lookups for better speed without changing the existing binning formulas
- signal is now reported as density per bp so regions with different bin widths are on the same scale

## Inputs

- `--co-bed`: BED file for CO signal
- `--dsb-bed`: BED file for DSB signal
- `--gff3`: gene annotation in GFF/GFF3 format
- `--gene-bed`: legacy BED-style gene annotation

Signal BED files must have at least 3 columns:
1. chromosome
2. start
3. end
4. count (optional)

If the count column is missing, the script uses `1`.

## Output

The default output file is `dual_anchor_metagene.pdf`.

The plotted values are density-normalized:
- upstream, downstream, TSS, and TTS bins are divided by their fixed bin width
- scaled gene-body bins are divided by their effective per-gene bin width (`body_length / number_of_body_bins`)
