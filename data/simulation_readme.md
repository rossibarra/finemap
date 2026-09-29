# Simulated Example Region Sets

These BED files each contain 100,000 simulated regions with lengths sampled from the empirical length distribution of `jri_v5.bed`.

All coordinates are on B73 v5 chromosomes `chr1`-`chr10`.

Sampling rules:

- `example1.bed`: uniform along the genome. Region starts are sampled uniformly across all valid genomic positions for the chosen length.
- `example2.bed`: enriched in 5' gene regions. One anchor point is sampled from strand-aware 5' 5 kb flanks, excluding positions that also fall in 3' flanks or the internal gene-body set below.
- `example3.bed`: enriched in gene bodies. One anchor point is sampled from internal gene-body segments after removing the first and last 5 kb of each gene; this excludes both 5' and 3' 5 kb flanks.
- `example4.bed`: enriched in 3' gene regions. One anchor point is sampled from strand-aware 3' 5 kb flanks, excluding positions that also fall in 5' flanks or the internal gene-body set.

Notes:

- Enrichment is defined by first drawing one point from the target annotation class, then assigning that point a random uniform position within the simulated region.
- Gene annotations come from `data/v5.genes.gff3`.
- Chromosome lengths come from `data/v5.fa.gz.fai`.
