# Make Chainfile: 

## Clone With Submodule

Clone this repository and initialize `AnchorWave/` in one step:

```bash
git clone --recurse-submodules <repo-url>
```

If you already cloned without submodules:

```bash
git submodule update --init --recursive
```

## Build AnchorWave (submodule)

Build the `AnchorWave` binary from the included submodule:

```bash
cd AnchorWave
cmake ./
make -j
cd ..
```

The alignment script will use `anchorwave` from your `PATH` if available, or fall back to the local binary under `AnchorWave/`.

## Files

v2 and v5 maize genome files, from [maizegdb](https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/):

- `v2.fa.gz`
- `v5.fa.gz`

You also need a **reference annotation file** (`.gff3` or `.gtf`) matching the reference genome used with `--ref`.

- `v5.gff.gz`

Scripts:

- `scripts/align_with_anchorwave.sh`

## Lift-Over to B73 v5

To move marker coordinates from AGPv2 to B73 v5, we converted the resulting MAF alignment from Anchorwave into a UCSC chain file with [maf-convert](https://gitlab.com/mcfrith/last/-/blob/main/doc/maf-convert.rst):

```bash
maf-convert chain AW_out/alignment.maf > v2v5.chain
```

## Requirements

Run `module load last` and `module load minimap2` (or install both tools so they are available on `PATH`).

## Quick Start

Run with explicit inputs:

```bash
scripts/align_with_anchorwave.sh \
  --gff path/to/reference.gff3 \
  --ref B73_RefGen_v2.fa.gz \
  --query Zm-B73-REFERENCE-NAM-5.0.fa.gz \
  --threads 16 \
  --outdir anchorwave_out
```

Typical run from repository root after building AnchorWave:

```bash
scripts/align_with_anchorwave.sh \
  --gff v5.gff3.gz \
  --ref v5.fa.gz \
  --query v2.fa.gz \
  --threads 16 \
  --outdir anchorwave_out
```

## Output

By default (`--outdir anchorwave_out`):

- `anchorwave_out/alignment.maf`
- `anchorwave_out/alignment.fragments.maf`
- `anchorwave_out/anchors.txt`
- `anchorwave_out/anchorwave.log`

## Notes for Maize-Scale Genomes

- Start with `--threads 16` (or higher based on available CPU).
- Use SSD-backed storage for faster SAM/MAF I/O.
- Keep sufficient free disk space: intermediate and final alignment files can be large.
- Default mode is `proali`; switch to `genoali` with `--mode genoali` if needed.

## Script Help

```bash
scripts/align_with_anchorwave.sh --help
```

# Ogut Map

## OGUT Marker Map Provenance

The OGUT marker map in this repository traces back to the marker resource reported in:

- Ogut F, Bian Y, Bradbury PJ, Holland JB. 2015. *Joint-multiple family linkage analysis predicts within-family variation better than single-family analysis of the maize nested association mapping population*. Heredity. [https://doi.org/10.1038/hdy.2014.123](https://doi.org/10.1038/hdy.2014.123)

The original marker table was recorded here as `ogut_fifthcM_map_agpv2.csv`, then converted to BED format as `ogut_fifthcM_map_agpv2.bed`. The BED conversion uses chromosome, physical position, marker name, source SNP ID, and cM value, with 0-based half-open coordinates:

```bash
awk -F',' 'NR>1 { pos=$4+0; print $3 "\t" pos-1 "\t" pos "\t" $2 "\t" $1 "\t" $5 }' \
  ogut_fifthcM_map_agpv2.csv > ogut_fifthcM_map_agpv2.bed
```

## Move to v5
After building the chain file, the marker BED was translated to v5 coordinates with [CrossMap](https://github.com/liguowang/CrossMap):

```bash
CrossMap bed v2v5.chain ogut_fifthcM_map_agpv2.bed ogut_v5.bed
```

This produced `ogut_v5.bed`, which contains the lifted marker positions on the v5 assembly while preserving the marker IDs and cM annotations. 

## Post-Lift Manual Marker Adjustments

After inspecting `ogut_v5.bed`, we found markers whose physical coordinates were not monotonically increasing within chromosome even though the cM values remained monotonic. To keep the map order consistent with the physical order on v5, we adjusted the cM values for the affected markers without changing their lifted physical positions.

Applied adjustments:

- On `Chr3`, `M2179` and `M2180` kept their physical coordinates, but their cM values were swapped.
- On `Chr7`, the block from `M5158` through `M5171` kept its physical coordinates, and the cM values across that block were reversed so that genetic order matches the physical order.

The pre-edit `ogut_v5.bed` was committed as `ce1ab1f`, and a local backup copy was also written to `ogut_v5.bed.bak` before those edits.

The final map is visualized in `ogut_v5_cm_vs_bp_by_chr.png`.

# Rodgers-Melnick Map

## Rodgers-Melnick Crossover Map Provenance

The Rodgers-Melnick crossover intervals in this repository trace back to:

- Rodgers-Melnick E, Bradbury PJ, Elshire RJ, Glaubitz JC, Acharya CB, Mitchell SE, Li C, Li Y, Buckler ES. 2015. *Recombination in diverse maize is stable, predictable, and associated with genetic load*. PNAS. [https://doi.org/10.1073/pnas.141386411](https://doi.org/10.1073/pnas.141386411)

The original interval tables were downloaded from Panzea and recorded here as:

- `RodgersMelnick2015PNAS_cnnamImputedXOsegments.txt`
- `RodgersMelnick2015PNAS_usnamImputedXOsegments.txt`

## Build the v2 BED with Unique IDs

We first removed heterozygous intervals, extracted chromosome/start/end plus family, and appended a unique ID to each v2 interval:

```bash
cat RodgersMelnick2015PNAS_cnnamImputedXOsegments.txt RodgersMelnick2015PNAS_usnamImputedXOsegments.txt | \
  grep -v 'het' | grep -v 'Family' | \
  tee >(cut -f 1,5,6 > left.txt) | cut -f 3 > right.txt

paste left.txt right.txt | \
  sed -e 's/\r//g' | \
  awk 'BEGIN{OFS="\t"} {print $0, sprintf("RMv2_%06d", NR)}' \
  > RodgersMelnickv2.bed
```

This produced `RodgersMelnickv2.bed` with columns:

- chromosome
- start
- end
- individual
- unique interval ID

## Add in Penny data

```cat penny_co_v2_male.txt penny_co_v2_female.txt | \
  grep -v '^chr[[:space:]]' | \
  sed -e 's/\r//g' | \
  awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $7, sprintf("Pennyv2_%06d", NR)}' \
  > Pennyv2.bed
  ```

## Add European maize data
```awk 'BEGIN{FS=OFS="\t"}
NR > 1 {
  print $3, $6, $7, $1, sprintf("EUROv2_%06d", NR-1)
}' hmm_co_events_long.tsv > eurov2.bed```

## Separate Marker Positions
```cat RodgersMelnickv2.bed eurov2.bed Pennyv2.bed | awk 'BEGIN{OFS="\t"} {
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


## Lift to v5

After building the chain file, the v2 interval BED was translated to B73 v5 coordinates with [CrossMap](https://github.com/liguowang/CrossMap):

```bash
CrossMap bed v2v5.chain markersv2.bed markersv5.bed
```

then rebuild into intervals
eep only column-5 IDs seen exactly twice total
drop IDs seen once or more than twice
drop IDs whose two hits are on different chromosomes
for kept IDs, output one merged BED interval using min start and max end
```awk '
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
' OFS='\t' markersv5.bed | sort -k1,1 -k2,2n >jri_v5.bed
```

Make 100kb bedfile
```awk 'BEGIN{FS=OFS="\t"; tile=100000}
$1 ~ /^chr([1-9]|10)$/ {
  chr = "Chr" substr($1, 4)
  len = $2
  for (start = 0; start < len; start += tile) {
    end = start + tile
    if (end > len) end = len
    print chr, start, end
  }
}' v5.fa.gz.fai > v5_100kb_tiles.bed```





## Convert Unique Regions into a Summed Rate Track

Each interval in `RMv5_unique.bed` represents one crossover localized to a region. To build a per-base crossover-rate track, we assigned each interval a uniform weight of `1 / interval_length`, emitted start and stop events, sorted them, and then swept across each chromosome to sum the active rates:

```bash
awk 'BEGIN{OFS="\t"}
{
  len = $3 - $2 + 1
  rate = 1 / len
  print $1, $2, rate
  print $1, $3 + 1, -rate
}' RMv5_unique.bed > RMv5.events.tsv

sort -k1,1V -k2,2n RMv5.events.tsv > RMv5.events.sorted.tsv

awk 'BEGIN{OFS="\t"; eps=1e-12}
function clean(x) {
  return (x < eps && x > -eps) ? 0 : x
}
NR==1 {
  chr = $1
  pos = $2
  delta = $3 + 0
  cur = 0
  next
}
$1 == chr && $2 == pos {
  delta += $3
  next
}
{
  cur = clean(cur + delta)

  if ($1 == chr && cur != 0 && pos < $2) {
    print chr, pos, $2 - 1, cur
  }

  if ($1 != chr) {
    cur = 0
  }

  chr = $1
  pos = $2
  delta = $3 + 0
}
END{
  cur = clean(cur + delta)
}' RMv5.events.sorted.tsv > RMv5_rates.bed
```

This sweep collapses all event deltas at the same position, applies the total change once, and then emits the constant-rate interval to the next event position. The final file, `RMv5_rates.bed`, is a chromosome-wise interval track in which each region stores the summed crossover rate over all uniquely resolved Rodgers-Melnick intervals lifted to B73 v5.

## Summarize Region CO and Chromosome cM Lengths

We next converted each `RMv5_rates.bed` interval into its total crossover contribution across the interval:

```bash
awk 'BEGIN{FS=OFS="\t"} {print $0, $4 * ($3 - $2 + 1)}' RMv5_rates.bed > RMv5_CO.bed
```

To anchor the Rodgers-Melnick map to the chromosome-scale genetic lengths in the Ogut v5 map, we computed chromosome cM spans from `ogut_v5.bed`:

```bash
awk 'BEGIN{FS=OFS="\t"}
{
  cm = $6 + 0
  if (!($1 in min) || cm < min[$1]) min[$1] = cm
  if (!($1 in max) || cm > max[$1]) max[$1] = cm
  n[$1]++
}
END{
  print "chromosome","cM_length","min_cM","max_cM","n_markers"
  for (chr in n) {
    print chr, max[chr] - min[chr], min[chr], max[chr], n[chr]
  }
}' ogut_v5.bed | sort -k1,1V > ogut_v5_chromosome_cM_lengths.tsv
```

## Rescale Rodgers-Melnick Regions to Ogut cM Lengths

Finally, we rescaled the Rodgers-Melnick interval contributions so that the summed crossover signal on each chromosome matches the total chromosome cM length from `ogut_v5_chromosome_cM_lengths.tsv`. We also computed `cM/Mb` for each region:

```bash
awk 'BEGIN{FS=OFS="\t"}
NR==FNR && FNR>1 { chr_cm[$1] = $2; next }
{
  len = $3 - $2 + 1
  val = $4 * len
  chr_sum[$1] += val
  row[NR] = $0
  row_val[NR] = val
  row_len[NR] = len
}
END{
  print "chr","start","end","rate","scaled_cM","cM_per_Mb"
  for (i = 1; i <= NR; i++) {
    split(row[i], f, FS)
    scaled_cm = (chr_sum[f[1]] > 0 ? row_val[i] * chr_cm[f[1]] / chr_sum[f[1]] : 0)
    cm_per_mb = (row_len[i] > 0 ? scaled_cm / row_len[i] * 1e6 : 0)
    print f[1], f[2], f[3], f[4], scaled_cm, cm_per_mb
  }
}' ogut_v5_chromosome_cM_lengths.tsv RMv5_rates.bed > RMv5.bed
```

This final `RMv5.bed` contains:

- chromosome
- start
- end
- summed Rodgers-Melnick crossover rate
- chromosome-scaled cM assigned to the interval
- `cM/Mb`

Comparison of the Rodgers-Melnick and Ogut maps can be seen in `RMv5_vs_ogut_cumulative_by_chromosome.png` (scaled so cumulative value is the same).
## Dual-Anchor Metagene Plotting

The script `metagene_dsb_co.py` plots one or two BED signal tracks around genes using two non-contiguous panels with explicit exon anchors:
- left panel: upstream flank, internal TSS window, first exon, last exon
- right panel: last exon, internal TTS window, downstream flank

The plot includes a visible axis break so it is clear that the full internal gene span between first and last exon is not shown continuously.

Current CLI shape:

```bash
python metagene_dsb_co.py   --gff v5.gff3   --bed1 example1.bed   --flank 10000 20   --gene 100 20
```

Arguments:
- `--gff` or `--gene_bed`: gene annotation input
- `--bed1`: first signal BED
- `--bed2`: optional second signal BED
- `--flank SIZE_BP N_BINS`: upstream and downstream flank size and bin count
- `--gene SIZE_BP N_BINS`: internal TSS/TTS boundary window size and bin count
- `--uniform`: spread each BED interval uniformly across its span instead of collapsing it to a midpoint

Current behavior:
- the full gene body is not plotted as a continuous scaled region
- first exon and last exon are shown as separate internal bins
- exon bins are normalized by exon length before averaging across genes
- if a gene has only one exon, that exon is used for both the first-exon and last-exon bins
- genes with length `<= 2 * INNER_SIZE` are excluded
- the script reports the excluded count to `stderr` as `X genes excluded because they are shorter than 2*INNER_SIZE`
- signal is normalized by per-gene bin width before averaging across genes, so unequal bin widths do not artificially depress narrow regions
- output is written to `dual_anchor_metagene.pdf`

