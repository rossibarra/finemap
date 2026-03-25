# Crossover Extraction Notes

## Input file

Source workbook:

- `/Users/jeffreyross-ibarra/Downloads/gb-2013-14-9-r103-S4.xlsx`

Sheet used:

- `Table_S3`

Columns used from the sheet:

- `map`
- `SNP_name`
- `chr_phy`
- `coordinate`
- `raw_data`

## Interpretation of the data

The `raw_data` field was treated as a genotype string for all individuals in a given mapping population (`map`) at a single SNP row.

Each character in `raw_data` was interpreted as one individual's genotype state at that SNP.

The workbook mixes multiple mapping populations, and the `raw_data` string length varies across `map` values. Because of that, crossover calling was done separately within each `map`.

The workbook does not expose individual names in this sheet, so anonymous sample IDs were assigned in the form:

- `MAP_ind001`
- `MAP_ind002`
- etc.

Examples:

- `CFD01_ind001`
- `CFF02_ind088`

## Base preprocessing steps

1. Read the `Table_S3` sheet from the Excel workbook.
2. Group all SNP rows by `map`.
3. Exclude rows with missing `map`, missing `SNP_name`, missing `raw_data`, missing physical coordinate, or `chr_phy == 0`.
4. Exclude short `raw_data` strings of length 6, which appear to be non-sample special cases rather than full genotype vectors.
5. Within each `map`, sort SNP rows by:
   - `chr_phy`
   - `coordinate`
   - `SNP_name`
6. For each character position in the `raw_data` string, treat that position as one anonymous individual.
7. Scan along each chromosome separately for each individual.
8. Ignore any genotype state that is not `A` or `B`.

## Unfiltered CO calling

The first-pass output (`results/co_events_long.tsv`) is the raw state-switch table.

In that file, a CO was called whenever two consecutive informative markers on the same chromosome changed from:

- `A` to `B`, or
- `B` to `A`

This means:

- `A -> B` was called as one CO.
- `B -> A` was called as one CO.
- `A -> A` was not called.
- `B -> B` was not called.
- Transitions through noninformative states such as `-`, `?`, `N`, `M`, `E`, or `#` were ignored until the next informative `A` or `B` was found.

## Filtered CO calling

To reduce likely false positives, a second-pass filtered caller was applied to produce:

- `results/co_events_long_filtered.tsv`

This filtered caller uses run-based haplotype support instead of accepting every single adjacent state change.

### Run definition

For each individual and chromosome, informative markers (`A` or `B`) were compressed into consecutive runs of the same state:

- `AAAA` becomes one `A` run
- `BBBB` becomes one `B` run
- a switch is evaluated only at a boundary between adjacent runs

Noninformative characters were skipped and did not contribute to run length.

### Filter rules used for a CO call

A filtered CO was called only when all of the following were true:

- The two flanking runs belonged to the same `map`.
- The two flanking runs were on the same physical chromosome (`chr_phy`).
- The earlier informative run and later informative run had different states.
- The new-state run contained at least `2` informative markers.
- The new-state run spanned at least `1,000,000` bp from its first to last informative marker.

In other words, a single marker flip was not enough. The new haplotype had to persist across multiple informative markers and across a nontrivial physical distance.

### Double-crossover cleanup

After the run-based calls were generated for each individual and chromosome, adjacent called COs were examined.

If two consecutive called COs for the same individual on the same chromosome were within:

- `2,000,000` bp

then both were removed as a likely short-interval double crossover artifact.

This specifically targets patterns like:

- `A -> B -> A` over a very short physical interval
- `B -> A -> B` over a very short physical interval

which are often caused by isolated genotype noise rather than true recombination.

## Event position definition

For both unfiltered and filtered calls, the CO position was reported as the midpoint between the physical coordinates of the two flanking markers at the accepted transition:

- left marker coordinate = `x1`
- right marker coordinate = `x2`
- reported CO position = `floor((x1 + x2) / 2)`

This is an inferred breakpoint interval midpoint, not an exact biological breakpoint.

## Output files

The workflow currently generates:

- `results/co_events_long.tsv`
- `results/co_events_long_filtered.tsv`
- `results/co_events_wide.tsv`

### `results/co_events_long.tsv`

One row per raw, unfiltered event, with:

- `sample_id`
- `map`
- `chromosome`
- `event_index`
- `event_bp`
- `left_coordinate`
- `right_coordinate`
- `left_snp`
- `right_snp`
- `left_genotype`
- `right_genotype`

### `results/co_events_long_filtered.tsv`

One row per filtered event, with the same columns as the unfiltered table:

- `sample_id`
- `map`
- `chromosome`
- `event_index`
- `event_bp`
- `left_coordinate`
- `right_coordinate`
- `left_snp`
- `right_snp`
- `left_genotype`
- `right_genotype`

### `results/co_events_wide.tsv`

Individuals are columns and rows contain successive event positions in the form:

- `chrN:position`

This wide table is still based on the raw, unfiltered event calls.

## Important caveats

- These are candidate crossover events based on genotype-state transitions, not validated biological breakpoints.
- The filtered table is more conservative than the raw table, but it is still heuristic.
- Marker-density differences among populations can affect how often runs meet the minimum-marker and minimum-span thresholds.
- The event position is bounded by the two flanking SNPs and should be interpreted as an approximate interval midpoint.
- The workbook still does not provide explicit individual names in this sheet, so all sample identifiers remain anonymous.

## Script used

The extraction was performed with:

- [scripts/extract_co_events.py](/Users/jeffreyross-ibarra/src/jri-arg/scripts/extract_co_events.py)

## Per-Individual Tiled Genome BED

A large per-individual BED file was also generated by combining:

- `jri_v5.bed`
- `v5_100kb_tiles.bed`

### Inputs used

#### `jri_v5.bed`

This file contains lifted crossover intervals with:

- chromosome
- start
- end
- individual identity
- interval ID

#### `v5_100kb_tiles.bed`

This file contains 100 kb genome tiles across the 10 maize v5 chromosomes:

- `Chr1` through `Chr10`

### Goal of the combined BED

For each individual in column 4 of `jri_v5.bed`, construct a genome-wide interval track that contains:

- the individual's original `jri_v5.bed` intervals
- all 100 kb tiles that do not overlap any of that individual's `jri_v5.bed` intervals

The result is a nonoverlapping BED for each individual, and then all individual-specific BEDs are concatenated into one final file.

### Overlap rule

For a given individual and chromosome:

- each original `jri_v5.bed` interval was kept as-is
- each 100 kb tile from `v5_100kb_tiles.bed` was tested against that individual's `jri_v5.bed` intervals
- if a tile overlapped any `jri_v5.bed` interval for that individual, the tile was dropped
- if a tile did not overlap any `jri_v5.bed` interval for that individual, the tile was retained

This means the final track for each individual contains no overlapping intervals between:

- original `jri_v5.bed` intervals
- retained background 100 kb tiles

### Output file

The concatenated output file is:

- `results/jri_v5_with_tiles.bed`

### Output columns

The output BED has 5 columns:

- column 1: chromosome
- column 2: start
- column 3: end
- column 4: individual identity
- column 5: interval class flag

Interpretation of column 5:

- `1` = original interval from `jri_v5.bed`
- `0` = retained nonoverlapping 100 kb tile from `v5_100kb_tiles.bed`

### Script used

This combined BED was generated with:

- `scripts/build_individual_tiled_track.py`

### Size of the final file

The final concatenated file is large because it expands the tiling across every individual represented in `jri_v5.bed`.

For the current run:

- output file: `results/jri_v5_with_tiles.bed`
- size: `9,762,317,527` bytes

This file is intended as a per-individual genome annotation track rather than a compact summary table.

## Gene Flank BED Generation

A combined gene-flank BED file was generated from:

- `v5.gff3`
- `data/v5.fa.gz.fai`

### Goal

Construct one combined BED file containing:

- all 2 kb regions 5' upstream of genes
- all 2 kb regions 3' downstream of genes

with strand-aware interpretation of upstream and downstream.

### Coordinate rules

Gene models were taken from `v5.gff3` using `gene` features only.

Because the GFF uses 1-based closed coordinates, intervals were converted to BED-style 0-based half-open coordinates before writing the output.

Chromosome lengths from `data/v5.fa.gz.fai` were used to clip flank intervals so they did not extend below 0 or beyond the chromosome end.

### Strand-aware definitions

For a `+` strand gene:

- 5' upstream = 2 kb immediately before the gene start
- 3' downstream = 2 kb immediately after the gene end

For a `-` strand gene:

- 5' upstream = 2 kb immediately after the gene end
- 3' downstream = 2 kb immediately before the gene start

### Output file

The combined flank BED file is:

- `results/v5_gene_flanks_2kb.bed`

### Output columns

This BED has 6 columns:

- column 1: chromosome
- column 2: start
- column 3: end
- column 4: gene ID
- column 5: flank class
- column 6: strand

Interpretation of column 5:

- `5prime_upstream`
- `3prime_downstream`

### Script used

This file was generated with:

- `scripts/extract_gene_flanks.py`

## Counting Gene-Flank Overlaps for the Per-Individual Tiled BED

To count, for each interval in `results/jri_v5_with_tiles.bed`, how many gene-flank intervals it overlaps, a two-step `bedtools intersect -c` pipeline was used.

### Chromosome-name normalization

`results/jri_v5_with_tiles.bed` uses chromosome names like:

- `Chr1`

while `results/v5_gene_flanks_2kb.bed` uses chromosome names like:

- `chr1`

So the tiled BED was normalized on the fly before intersection:

- `Chr1` -> `chr1`
- `Chr2` -> `chr2`
- etc.

### Bedtools command

```bash
awk 'BEGIN{FS=OFS="\t"}
{
  chrom = $1
  if (chrom ~ /^Chr/) chrom = "chr" substr(chrom, 4)
  print chrom, $2, $3, $4, $5
}' results/jri_v5_with_tiles.bed \
| bedtools intersect -a stdin \
    -b <(awk '$5=="5prime_upstream"' results/v5_gene_flanks_2kb.bed) \
    -c \
| bedtools intersect -a stdin \
    -b <(awk '$5=="3prime_downstream"' results/v5_gene_flanks_2kb.bed) \
    -c \
> results/jri_v5_with_tiles_flank_counts.bed
```

### Resulting output columns

The resulting file has 7 columns:

- column 1: chromosome
- column 2: start
- column 3: end
- column 4: individual identity
- column 5: original `0/1` interval flag
- column 6: count of overlapping `5prime_upstream` intervals
- column 7: count of overlapping `3prime_downstream` intervals
