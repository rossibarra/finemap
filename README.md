# finemap

This repository contains workflows for deriving maize crossover interval sets, lifting them onto B73 v5 coordinates, generating simulated comparison region sets, and plotting gene-centered summary tracks.

The material below consolidates the previous `euro_readme.md`, `meta_readme.md`, `simulation_readme.md`, and `hmm_methods.md` notes into one project README, with the workflow split into:

- data manipulation and lift-over
- simulation data
- plotting code

## Data Manipulation And Lift-Over

### European Workbook Crossover Extraction

European crossover candidates were extracted from the workbook:

- `/Users/jeffreyross-ibarra/Downloads/gb-2013-14-9-r103-S4.xlsx`

Sheet used:

- `Table_S3`

Columns used:

- `map`
- `SNP_name`
- `chr_phy`
- `coordinate`
- `raw_data`

`raw_data` was interpreted as a per-population genotype string, where each character corresponds to one anonymous individual at one SNP. Because the workbook mixes populations and the string length varies by `map`, all calling was done separately within each population.

Anonymous sample IDs were assigned as:

- `MAP_ind001`
- `MAP_ind002`
- etc.

Base preprocessing:

1. Read `Table_S3`.
2. Group rows by `map`.
3. Drop rows with missing `map`, missing `SNP_name`, missing `raw_data`, missing coordinate, or `chr_phy == 0`.
4. Drop short `raw_data` strings of length `6`.
5. Sort each population by chromosome and coordinate.
6. Treat each character position in `raw_data` as one individual.
7. Scan each chromosome independently.
8. Ignore genotype states other than `A` and `B`.

#### Raw Switch-Based Calls

The raw event table is:

- `results/co_events_long.tsv`

A crossover was called whenever consecutive informative markers changed from `A -> B` or `B -> A`. Noninformative states were skipped until the next informative marker.

#### Filtered Switch-Based Calls

The filtered event table is:

- `results/co_events_long_filtered.tsv`

Filtering was run on informative-state runs rather than single-marker flips. A filtered call required:

- a state change between adjacent informative runs
- at least `2` informative markers in the new-state run
- at least `1,000,000` bp spanned by the new-state run

After run-based calling, adjacent calls within `2,000,000` bp for the same sample and chromosome were removed as likely short double-crossover artifacts.

For both raw and filtered calls, `event_bp` is:

- `floor((left_coordinate + right_coordinate) / 2)`

Scripts used:

- `scripts/extract_co_events.py`

### HMM-Based Crossover Workflow

The preferred crossover-calling workflow is the HMM pipeline, which uses the same workbook and sheet as above but applies per-population marker cleaning and chromosome-wise decoding before interval calling.

Main outputs:

- `results/hmm_cleaned_matrices/`
- `results/hmm_co_events_long.tsv`
- `results/hmm_qc_summary.tsv`
- `results/marey_map_hmm_co_events.png`

#### Marker Cleaning

Cleaning was done separately within each `map`:

1. Remove rows with missing required fields, `chr_phy == 0`, or `raw_data` length `<= 6`.
2. Deduplicate markers sharing the same `(chromosome, coordinate)` by keeping the row with the highest informative genotype count.
3. Remove individuals with missingness greater than `0.30`.
4. Remove markers with missingness greater than `0.20`.
5. Remove markers with informative `A` fraction outside `0.10` to `0.90`.
6. Remove markers with isolated-flip rate greater than `0.12`.

This run retained:

- `23` cleaned population matrices
- `2,209` individuals

#### HMM Model

The HMM was run separately for each population, chromosome, and individual.

States:

- `A`
- `B`

Observations:

- `A`
- `B`
- missing/noninformative

Emission probabilities:

- match probability `0.98`
- mismatch probability `0.02`
- missing probability `0.50` under either state

Distance-aware transition settings:

- base rate per bp `1e-8`
- minimum transition probability `1e-6`
- maximum transition probability `0.05`

The hidden path was decoded with Viterbi.

#### HMM CO Interval Calling

A crossover interval was called whenever the decoded state changed between adjacent markers on the same chromosome. Each event records:

- `sample_id`
- `map`
- `chromosome`
- `event_index`
- `event_bp`
- `left_coordinate`
- `right_coordinate`
- `left_snp`
- `right_snp`
- `left_state`
- `right_state`

`event_bp` is again the midpoint:

- `floor((left_coordinate + right_coordinate) / 2)`

After HMM calling, adjacent events within `2,000,000` bp for the same individual and chromosome were removed.

Run summary for the recorded execution:

- `32,439` HMM-based CO intervals
- `23` cleaned population matrices
- `2,209` retained individuals with HMM calls

Scripts used:

- `scripts/hmm_co_pipeline.py`
- `scripts/plot_marey_map.py`

### Lift-Over To B73 v5

Several workflows in this repository rely on lifting marker or interval coordinates from AGPv2 to B73 v5.

#### Chain File Construction

The v2-to-v5 chain file was derived from an AnchorWave whole-genome alignment.

Build the submodule:

```bash
git clone --recurse-submodules <repo-url>
cd AnchorWave
cmake ./
make -j
cd ..
```

Run the alignment:

```bash
scripts/align_with_anchorwave.sh \
  --gff v5.gff3.gz \
  --ref v5.fa.gz \
  --query v2.fa.gz \
  --threads 16 \
  --outdir anchorwave_out
```

Convert the resulting MAF to a UCSC chain:

```bash
maf-convert chain anchorwave_out/alignment.maf > v2v5.chain
```

Requirements:

- `last`
- `minimap2`
- `AnchorWave`

Primary alignment script:

- `scripts/align_with_anchorwave.sh`

#### Rodgers-Melnick, Penny, And European Interval Inputs

The repository combines multiple v2-era crossover interval sources before lift-over:

- Rodgers-Melnick NAM intervals
- Penny interval data
- European workbook-derived events

Rodgers-Melnick source:

- Rodgers-Melnick E, Bradbury PJ, Elshire RJ, Glaubitz JC, Acharya CB, Mitchell SE, Li C, Li Y, Buckler ES. 2015. *Recombination in diverse maize is stable, predictable, and associated with genetic load*. PNAS.
- downloaded interval tables:
  `RodgersMelnick2015PNAS_cnnamImputedXOsegments.txt`
  and
  `RodgersMelnick2015PNAS_usnamImputedXOsegments.txt`
- described in the earlier repository notes as downloaded from Panzea

Penny source:

- Kianian PMA, Wang M, Simons K, et al. 2018. *High-resolution crossover mapping reveals similarities and differences of male and female recombination in maize*. Nature Communications. `https://doi.org/10.1038/s41467-018-04562-5`
- local interval tables:
  `penny_co_v2_male.txt`
  and
  `penny_co_v2_female.txt`

European source:

- Bauer E, Falque M, Walter H, et al. 2013. *Intraspecific variation of recombination rate in maize*. Genome Biology. `https://doi.org/10.1186/gb-2013-14-9-r103`
- workbook: `/Users/jeffreyross-ibarra/Downloads/gb-2013-14-9-r103-S4.xlsx`
- sheet: `Table_S3`
- columns used for crossover extraction:
  `map`, `SNP_name`, `chr_phy`, `coordinate`, and `raw_data`

Example conversions to v2 BED:

```bash
cat RodgersMelnick2015PNAS_cnnamImputedXOsegments.txt RodgersMelnick2015PNAS_usnamImputedXOsegments.txt | \
  grep -v 'het' | grep -v 'Family' | \
  tee >(cut -f 1,5,6 > left.txt) | cut -f 3 > right.txt

paste left.txt right.txt | \
  sed -e 's/\r//g' | \
  awk 'BEGIN{OFS="\t"} {print $0, sprintf("RMv2_%06d", NR)}' \
  > RodgersMelnickv2.bed
```

```bash
cat penny_co_v2_male.txt penny_co_v2_female.txt | \
  grep -v '^chr[[:space:]]' | \
  sed -e 's/\r//g' | \
  awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $7, sprintf("Pennyv2_%06d", NR)}' \
  > Pennyv2.bed
```

```bash
awk 'BEGIN{FS=OFS="\t"} NR > 1 {
  print $3, $6, $7, $1, sprintf("EUROv2_%06d", NR-1)
}' hmm_co_events_long.tsv > eurov2.bed
```

The interval endpoints were split into marker positions before lift-over:

```bash
cat RodgersMelnickv2.bed eurov2.bed Pennyv2.bed | awk 'BEGIN{OFS="\t"} {
  x = $2
  y = $3 - 1

  $2 = x
  $3 = x + 1
  print

  $2 = y
  $3 = y + 1
  print
}' > markersv2.bed
```

Lift to v5:

```bash
CrossMap bed v2v5.chain markersv2.bed markersv5.bed
```

#### `jri_v5.bed` Creation

`jri_v5.bed` is the combined lifted interval set derived from the Rodgers-Melnick, Penny, and European inputs above. It is not a raw concatenation of the original interval tables.

Construction steps:

1. Build `RodgersMelnickv2.bed`, `Pennyv2.bed`, and `eurov2.bed`.
2. Concatenate those v2 interval BEDs.
3. Split each interval into two single-base endpoint markers to create `markersv2.bed`.
4. Lift those endpoints to v5 with CrossMap, producing `markersv5.bed`.
5. Rebuild intervals from the lifted endpoint pairs.

During interval rebuilding:

- keep only IDs seen exactly twice total after lift-over
- drop IDs seen once or more than twice
- drop IDs whose two lifted endpoints land on different chromosomes
- for retained IDs, emit one interval using the minimum lifted start and maximum lifted end

The earlier repository notes recorded the interval-rebuild step as:

```bash
awk '
{
  id = $5
  chrkey = $1 FS id

  total[id]++
  perchr[chrkey]++

  if (!(chrkey in min_start) || $2 < min_start[chrkey]) min_start[chrkey] = $2
  if (!(chrkey in max_end)   || $3 > max_end[chrkey])   max_end[chrkey]   = $3
}
END {
  for (chrkey in perchr) {
    split(chrkey, a, FS)
    chr = a[1]
    id  = a[2]

    if (total[id] == 2 && perchr[chrkey] == 2) {
      print chr, min_start[chrkey], max_end[chrkey], indiv, id
    }
  }
}
' OFS='\t' markersv5.bed | sort -k1,1 -k2,2n > jri_v5.bed
```

This produces a v5 BED of lifted crossover intervals with one row per successfully reconstructed interval.

#### Ogut Marker Map

The Ogut map originates from:

- Ogut F, Bian Y, Bradbury PJ, Holland JB. 2015. *Joint-multiple family linkage analysis predicts within-family variation better than single-family analysis of the maize nested association mapping population*. Heredity.

Source files:

- `ogut_fifthcM_map_agpv2.csv`
- `ogut_fifthcM_map_agpv2.bed`

Lift-over:

```bash
CrossMap bed v2v5.chain ogut_fifthcM_map_agpv2.bed ogut_v5.bed
```

Manual post-lift cM adjustments were applied to preserve chromosome-wise monotonicity:

- on `Chr3`, `M2179` and `M2180` had their cM values swapped
- on `Chr7`, markers `M5158` through `M5171` had their cM values reversed across the block

#### Rodgers-Melnick Rate Track

After reconstructing unique v5 intervals, the Rodgers-Melnick intervals were converted into a per-base rate track by assigning each interval a weight of `1 / interval_length`, sweeping start and stop events, and then rescaling per chromosome to match Ogut chromosome cM spans.

Main derived outputs:

- `RMv5_rates.bed`
- `RMv5_CO.bed`
- `ogut_v5_chromosome_cM_lengths.tsv`
- `RMv5.bed`

## Simulation Data

The repository includes simulated BED region sets with lengths sampled from the empirical length distribution of `jri_v5.bed`.

Files:

- `example1.bed`
- `example2.bed`
- `example3.bed`
- `example4.bed`

All coordinates are on B73 v5 chromosomes `chr1` through `chr10`.

Sampling rules:

- `example1.bed`: uniform genome-wide placement
- `example2.bed`: enriched in strand-aware 5' gene flanks
- `example3.bed`: enriched in internal gene-body segments after removing the first and last 5 kb
- `example4.bed`: enriched in strand-aware 3' gene flanks

Notes:

- enrichment is defined by sampling an anchor point from the target annotation class and then placing the simulated interval uniformly around that point
- gene annotations come from `v5.gff3`
- chromosome lengths come from `data/v5.fa.gz.fai`

## Plotting Code

### `metaplot.py`

`metaplot.py` builds a gene-centered metaplot from a GFF/GFF3 annotation and an input BED file.

Behavior:

- reads `gene` features from `--gff`
- normalizes chromosome names
- truncates flanks so they do not extend into neighboring genes
- plots mean signal across windows
- adds a standard-error ribbon
- adds a lower panel showing contributing gene counts
- does not explicitly exclude genes by length
- reports chromosome-name normalization warnings to `stderr` when renaming occurs

Signal modes:

- default midpoint mode: each BED interval contributes at its midpoint
- `--uniform`: distribute interval value across its full span

BED input expectations:

1. chromosome
2. start
3. end
4. value optional

If column 4 is missing or nonnumeric, value defaults to `1`.

Basic usage:

```bash
python metaplot.py \
  --gff v5.gff3 \
  --input jri_small.bed \
  --bin-size 100 \
  --flanking-bp 1000 \
  --body-bins 25
```

```bash
python metaplot.py \
  --gff v5.gff3 \
  --input jri_small.bed \
  --bin-size 100 \
  --flanking-bp 1000 \
  --body-bins 25 \
  --uniform \
  --output my_metaplot.pdf \
  --no-show
```

Main arguments:

- `--gff`
- `--input`
- `--bin-size`
- `--flanking-bp`
- `--body-bins`
- `--uniform`
- `--output`
- `--title`
- `--no-show`

Default output:

- `metaplot.pdf`

Display behavior:

- the plot is saved as a PDF
- unless `--no-show` is used, the plot is also shown interactively

### Marey-Style Plotting

The HMM workflow also produces a Marey-style recombination summary plot with:

- `scripts/plot_marey_map.py`

Behavior:

- accepts optional input and output paths as positional arguments
- defaults to reading `results/co_events_long.tsv`
- defaults to writing `results/marey_map_co_events.png`
- writes a PNG figure with `fig.savefig(...)`
- prints `Wrote <output_path>` to standard output

## Notes

- The HMM-based workflow should be treated as the preferred source of final crossover intervals.
- The older switch-based European tables remain useful as an exploratory and provenance-preserving intermediate.
- Interval locations remain bounded by marker spacing, so all reported crossover positions are interval midpoints rather than exact breakpoints.
