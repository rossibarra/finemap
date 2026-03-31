# metaplot.py

`metaplot.py` builds a gene-centered metaplot from a GFF gene annotation and an input BED file.

## What it does

- reads regions labeled `gene` from `--gff`
- uses the midpoint of each interval in `--input`
- sums midpoint values into fixed-width bins
- treats missing or non-numeric BED column 4 as `1`
- truncates 5' and 3' flanks so they do not extend into neighboring gene bodies
- averages each plotted bin over only the genes that contribute data to that bin

Inside genes, the plot shows:
- up to `--body-bins` bins immediately 3' of the TSS
- up to `--body-bins` bins immediately 5' of the TTS

If a gene is shorter than `2 * --body-bins * --bin-size`, the script uses `floor(gene_length / (2 * bin_size))` bins on each side and drops any leftover central segment that does not fill a complete bin pair.

## Basic usage

```bash
python metaplot.py \
  --gff v5.gff3 \
  --input jri_small.bed \
  --bin-size 100 \
  --flanking-bp 1000 \
  --body-bins 25
```

## Main arguments

- `--gff`: GFF/GFF3 file containing `gene` features
- `--input`: BED file with 3 or 4 columns
- `--bin-size`: bin width in bp
- `--flanking-bp`: maximum 5' and 3' flank size in bp
- `--body-bins`: number of internal bins plotted on each side of the gene
- `--output`: output PDF path, default `metaplot.pdf`
- `--no-show`: save the plot without opening a window
