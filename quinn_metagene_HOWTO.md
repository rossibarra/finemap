# quinn_metagene.py how-to

## Typical run with GFF3

```bash
python quinn_metagene.py \
  --gff3 data/v5.genes.gff3 \
  --co-bed co.bed \
  --dsb-bed dsb.bed \
  --flank-bin-size 500 \
  --body-bins 50 \
  --inner-bin-size 100
```

## Typical run with the legacy gene BED

```bash
python quinn_metagene.py \
  --gene-bed genes.bed \
  --co-bed co.bed \
  --dsb-bed dsb.bed \
  --flank-bin-size 500 \
  --body-bins 50 \
  --inner-bin-size 100
```

## Useful optional flags

- `--flank-size 10000`: total upstream and downstream flank size in bp
- `--inner-size 2000`: TSS and TTS inner-window size in bp
- `--smooth-window 5`: moving-average smoothing window
- `--output custom.pdf`: output PDF name
- `--no-show`: skip the interactive plot window

## Notes

- chromosome names are normalized so `Chr1` and `chr1` match
- GFF3 input uses `gene` features and converts starts to 0-based coordinates before plotting
- if your BED files have only 3 columns, each interval contributes a count of `1`
- output values are density per bp, so differing flank, inner, and scaled body bin widths are normalized before averaging
