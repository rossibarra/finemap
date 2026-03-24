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
