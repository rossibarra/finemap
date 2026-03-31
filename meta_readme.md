# metaplot.py

`metaplot.py` builds a gene-centered metaplot from a GFF/GFF3 gene annotation and an input BED file.

## What it does

- reads features labeled `gene` from `--gff`
- normalizes chromosome names so `Chr1` and `chr1` match
- reports chromosome-name normalizations to `stderr`
- truncates 5' and 3' flanks so they do not extend into neighboring gene bodies
- plots the average signal per window across genes
- draws a semi-transparent standard-error band around the mean line
- draws a lower histogram panel showing how many genes contribute to each plotted window

## Signal modes

By default, `metaplot.py` uses midpoint mode:

- each BED interval is reduced to its midpoint
- the value of that interval is added to the single plot window containing the midpoint

With `--uniform`:

- each BED interval is treated as uniformly distributed across its full span
- the interval contributes `value / interval_length` per bp
- each plot window receives the overlap-weighted sum from all intervals overlapping that window

## Input BED rules

The BED file given with `--input` must have at least 3 columns:

1. chromosome
2. start
3. end
4. value (optional)

If column 4 is missing, each interval is assigned a value of `1`.

If column 4 is present but non-numeric, that row is also assigned a value of `1`.

## Plot layout

The plot is built from four pieces:

- upstream flank up to `--flanking-bp`
- `--body-bins` windows immediately 3' of the TSS
- a shaded blank gap separating the TSS-side and TTS-side gene windows
- `--body-bins` windows immediately 5' of the TTS
- downstream flank up to `--flanking-bp`

Each gene contributes only the windows that exist after flank truncation and gene-length limits are applied.

For short genes:

- if the gene is shorter than `2 * --body-bins * --bin-size`
- the script uses `floor(gene_length / (2 * bin_size))` bins on each side
- those TSS-side bins are anchored to the TSS
- those TTS-side bins are anchored to the TTS
- any leftover center segment that does not fill a complete bin pair is dropped

Bins that are truncated away for a given gene do not contribute to the mean, standard error, or lower count panel for that plot position.

## Basic usage

Midpoint mode:

```bash
python metaplot.py \
  --gff v5.gff3 \
  --input jri_small.bed \
  --bin-size 100 \
  --flanking-bp 1000 \
  --body-bins 25
```

Uniform mode:

```bash
python metaplot.py \
  --gff v5.gff3 \
  --input jri_small.bed \
  --bin-size 100 \
  --flanking-bp 1000 \
  --body-bins 25 \
  --uniform
```

Non-interactive run:

```bash
python metaplot.py \
  --gff v5.gff3 \
  --input jri_small.bed \
  --bin-size 100 \
  --flanking-bp 1000 \
  --body-bins 25 \
  --output my_metaplot.pdf \
  --no-show
```

## Main arguments

- `--gff`: GFF/GFF3 file containing `gene` features
- `--input`: BED file with 3 or 4 columns
- `--bin-size`: bin width in bp
- `--flanking-bp`: maximum 5' and 3' flank size in bp
- `--body-bins`: number of internal bins plotted on each side of the gene
- `--uniform`: distribute each interval across its full span instead of using only the midpoint
- `--output`: output PDF path, default `metaplot.pdf`
- `--title`: plot title, default `Metaplot`
- `--no-show`: save the plot without opening a window

## Output

The script writes a PDF plot, by default:

```text
metaplot.pdf
```

The figure contains:

- an upper panel with the mean signal and standard-error band
- a lower panel with the number of genes contributing to each window
