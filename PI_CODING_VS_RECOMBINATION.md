# Coding-Site Diversity (0D / 4D) vs Recombination Rate

The coding-sequence version of [PI_VS_RECOMBINATION.md](PI_VS_RECOMBINATION.md). Nucleotide
diversity (π) is measured separately at **0-fold degenerate** sites (a nonsynonymous proxy) and
**4-fold degenerate** sites (a synonymous proxy) in the same 100 kb windows as the FineMap
hierarchical map, in maize, *Zea mays* ssp. *mexicana*, and *Zea mays* ssp. *parviglumis*.

Splitting diversity this way adds something the all-sites analysis cannot give: 4D sites are
close to neutral, so π<sub>4D</sub> tracks the local effective population size, while 0D sites
are under purifying selection. Their ratio π<sub>0D</sub>/π<sub>4D</sub> is a within-window
measure of how efficiently deleterious variation is removed, and under linked selection it
should *fall* as recombination rises — selection is more efficient where linkage with
neighbouring selected sites is broken up faster.

## Table of Contents

- [Input Data](#input-data)
- [Why 0D/4D Rather Than Literal Synonymous/Nonsynonymous](#why-0d4d-rather-than-literal-synonymousnonsynonymous)
- [Degeneracy Annotation](#degeneracy-annotation)
- [Estimator and Windows](#estimator-and-windows)
- [Running the Analysis](#running-the-analysis)
  - [Step 1 — Reference Genome and Annotation](#step-1--reference-genome-and-annotation)
  - [Step 2 — 0D/4D Site Annotation](#step-2--0d4d-site-annotation)
  - [Step 3 — Windowed π at 0D and 4D Sites](#step-3--windowed--at-0d-and-4d-sites)
  - [Step 4 — Correlation and Figures](#step-4--correlation-and-figures)
- [Callable-Site Threshold](#callable-site-threshold)
- [Results](#results)
- [Caveats](#caveats)

## Input Data

The same per-chromosome all-sites VCFs as the all-sites analysis: 29 alignment-derived haploid
pseudo-genomes — 8 maize inbreds, 12 *mexicana*, 9 *parviglumis* — with invariant sites retained
so that callable-site counts are known.

**These VCFs are unpublished and are not distributed with this repository.** Neither are the π
tables derived from them; `data/pixy/` is git-ignored.

Two reference files are downloaded from MaizeGDB (both git-ignored, see
[Step 1](#step-1--reference-genome-and-annotation)):

| File | Size | Purpose |
|------|------|---------|
| `Zm-B73-REFERENCE-NAM-5.0.fa.gz` | 646 MB | reference bases, to reconstruct codons |
| `Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz` | 12 MB | CDS features, strand and phase |

The repository's own `data/v5.genes.gff3` contains only `gene` features and cannot be used here;
it is still used as the gene-density covariate, exactly as in the all-sites analysis. The
downloaded FASTA's sequence names (`chr1`…`chr10`) and lengths are byte-identical to
`data/v5.fa.gz.fai`, so no coordinate or naming translation is needed against the VCFs.

## Why 0D/4D Rather Than Literal Synonymous/Nonsynonymous

Calling a *variant* synonymous or nonsynonymous requires knowing which base changed and in what
codon context, and the answer can differ between the two alleles at a multiallelic site. Worse,
it makes the denominator ill-defined: an unbiased π needs the number of *callable* synonymous
positions, not just the number of synonymous SNPs, and a site is only partly synonymous when
two of its three possible changes are silent and one is not.

Degeneracy classes remove both problems:

- A **4-fold degenerate** position is one where all three alternative bases encode the same
  amino acid. Every mutation there is synonymous, so both numerator and denominator are clean.
- A **0-fold degenerate** position is one where all three alternatives change the amino acid.
  Every mutation there is nonsynonymous.

2-fold and 3-fold positions are discarded outright. The cost is throwing away roughly a third of
coding positions; the benefit is that π<sub>4D</sub> and π<sub>0D</sub> are per-site diversities
over well-defined, mutually exclusive site classes and are directly comparable to each other.

Caveat inherent to the proxy: degeneracy is assigned from the **B73 reference codon**. At a site
that is itself polymorphic the alternative allele may sit in a codon of a different degeneracy
class, so a small fraction of sites are misclassified relative to a full per-allele annotation.

## Degeneracy Annotation

`scripts/degenerate_sites.py` builds the annotation:

1. **One transcript per gene.** Every mRNA with `biotype=protein_coding` is read from the GFF3.
   Per gene the transcript flagged `canonical_transcript=1` is used, falling back to the longest
   CDS if none is flagged. In this annotation all 39,035 protein-coding genes have a flagged
   canonical transcript, chosen out of 71,791 protein-coding transcripts. Using one transcript
   per gene is what stops a position being counted twice through alternative isoforms.
2. **Codon reconstruction.** CDS features are concatenated in transcript order (descending
   genomic coordinate and reverse-complemented for `-` strand genes), the leading partial codon
   indicated by the first CDS's GFF phase is trimmed, the trailing remainder is dropped, and the
   rest is reshaped into codons. Codons containing any non-ACGT base are skipped, as are stop
   codons.
3. **Classification.** A 64 × 3 lookup table built from the standard genetic code assigns each
   codon position 4-fold, 0-fold, or neither, and the whole transcript is classified in one
   vectorised pass.
4. **Deduplication.** Positions claimed by more than one gene are collapsed; the 205,113
   positions given *conflicting* classes by overlapping genes (usually genes on opposite
   strands) are dropped from both classes.

Output is `data/degenerate_sites_v5.npz`, holding per chromosome two sorted `int32` arrays of
1-based positions, `<chrom>_0D` and `<chrom>_4D` (12 MB, git-ignored). The whole annotation
takes about 15 s.

**Counts and sanity check** (chr1–chr10):

| Class | Sites | Codon position 1 | 2 | 3 |
|-------|-------|------------------|---|---|
| 0D | 27,202,315 | 46.52% | 51.61% | 1.87% |
| 4D | 7,179,256 | 0.00% | 0.00% | 100.00% |

This is exactly the expected pattern. 4D sites are *only* third positions. 0D sites are
overwhelmingly first and second positions; the 1.87% at third positions are real — the third
positions of ATG (Met) and TGG (Trp), the two codons with no synonymous third-position change.
The 0D:4D ratio of 3.8 also matches the standard code, in which 0-fold positions outnumber
4-fold ones roughly four to one.

## Estimator and Windows

Identical to the all-sites analysis. Per site, with *n* non-missing haplotypes and *c<sub>a</sub>*
the count of allele *a*:

```
diffs_site = (n² - Σ c_a²) / 2
pairs_site = n (n - 1) / 2
```

summed separately over the sites in a window, with π the ratio of the sums. Invariant sites
contribute to the denominator only — which is why the all-sites VCF is essential; sites with
*n* < 2 contribute to neither. The only change is that a site must also appear in the 0D (or 4D)
position list to contribute at all: `scripts/haploid_pi.py` gained a `--restrict` option that
tests VCF positions against the sorted per-chromosome array with `np.searchsorted`, so no
per-chromosome mask is ever materialised and the whole-genome pass still takes a few minutes.

Windows are the 21,325 uniform 100 kb windows of `data/finemap_hierarchical_v5.bed`, and
`--lowercase-chrom` reconciles the BED's `Chr1` with the VCF's `chr1`.

Of the annotated sites, 22,002,658 of 27,202,315 0D positions (80.9%) and 5,802,551 of 7,179,256
4D positions (80.8%) are callable in at least two haplotypes somewhere in the VCFs.

## Running the Analysis

### Step 1 — Reference Genome and Annotation

```bash
cd data
B=https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0
curl -O $B/Zm-B73-REFERENCE-NAM-5.0.fa.gz
curl -O $B/Zm-B73-REFERENCE-NAM-5.0.fa.gz.fai
curl -O $B/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
# confirm the FASTA matches the VCF's chromosome naming and lengths
diff <(cut -f1,2 Zm-B73-REFERENCE-NAM-5.0.fa.gz.fai) <(cut -f1,2 v5.fa.gz.fai)
```

These files are large and are listed in `.gitignore`; MaizeGDB is preferred over Ensembl Plants
because Ensembl names the chromosomes `1`…`10` rather than `chr1`…`chr10`.

### Step 2 — 0D/4D Site Annotation

```bash
python scripts/degenerate_sites.py \
  --fasta data/Zm-B73-REFERENCE-NAM-5.0.fa.gz \
  --gff data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz \
  --out data/degenerate_sites_v5.npz
```

| Argument | Description |
|----------|-------------|
| `--fasta` | Reference FASTA, gzipped or plain |
| `--gff` | Full GFF3 **with CDS features**; `data/v5.genes.gff3` will not work |
| `--chroms` | Sequences to annotate, default `chr1`…`chr10` |
| `--out` | Output `.npz` of `<chrom>_0D` / `<chrom>_4D` position arrays |

### Step 3 — Windowed π at 0D and 4D Sites

```bash
for CLS in 0D 4D; do
  python scripts/haploid_pi.py \
    --vcf combined.chr{1..10}.all_sites.vcf.gz \
    --populations data/pixy/populations.txt \
    --windows data/finemap_hierarchical_v5.bed \
    --lowercase-chrom \
    --restrict data/degenerate_sites_v5.npz --restrict-key $CLS \
    --out data/pixy/pi_${CLS}_finemap_windows.tsv
done
```

| Argument | Description |
|----------|-------------|
| `--restrict` | `.npz` of `<chrom>_<key>` 1-based position arrays; only those sites count |
| `--restrict-key` | Array suffix to use, `0D` or `4D` (default `4D`) |

All other arguments are as documented in [PI_VS_RECOMBINATION.md](PI_VS_RECOMBINATION.md).
Output columns are pixy's: `pop`, `chromosome`, `window_pos_1`, `window_pos_2`, `avg_pi`,
`no_sites`, `count_diffs`, `count_comparisons`.

### Step 4 — Correlation and Figures

```bash
python scripts/pi_coding_vs_recombination.py \
  --pi-0d data/pixy/pi_0D_finemap_windows.tsv \
  --pi-4d data/pixy/pi_4D_finemap_windows.tsv \
  --bed data/finemap_hierarchical_v5.bed \
  --gff data/v5.genes.gff3 \
  --centromeres data/NAM_centromere_coords-cenH3.csv \
  --min-sites-0d 200 --min-sites-4d 100 \
  --out-prefix results/pi_coding_vs_recombination
```

| Argument | Description |
|----------|-------------|
| `--pi-0d`, `--pi-4d` | pixy-format π tables from Step 3 |
| `--bed` | FineMap BED; column 6 (cM/Mb) is the rate |
| `--gff`, `--centromeres` | Covariates for the partial correlation |
| `--min-sites-0d` | Minimum callable 0D sites per window, default 200 |
| `--min-sites-4d` | Minimum callable 4D sites per window, default 100 |
| `--block-mb` | Bootstrap block size in Mb, default 5 |
| `--n-boot` | Bootstrap replicates, default 1000 |

Spearman ρ, the moving-block bootstrap and the partial correlation are imported directly from
`scripts/pi_vs_recombination.py`, so the two analyses use identical statistics. Outputs are
`results/pi_coding_vs_recombination.png`, `_ratio.png`, `_summary.tsv`, `_deciles.tsv` and
`_windows.tsv.gz`.

## Callable-Site Threshold

The all-sites analysis used `--min-sites 5000`, which is impossible here: coding sequence is a
few percent of the genome, and the median 100 kb window holds only 594 callable 0D sites and 159
callable 4D sites (a quarter of windows hold none at all, having no genes). Thresholds were set
to keep essentially every window that contains genes at all, while excluding windows whose π
would rest on a handful of sites:

| Threshold | 0D windows kept | Threshold | 4D windows kept |
|-----------|-----------------|-----------|-----------------|
| ≥ 20 | 13,195 (61.9%) | ≥ 20 | 13,144 (61.6%) |
| ≥ 50 | 13,173 (61.8%) | ≥ 50 | 12,879 (60.4%) |
| ≥ 100 | 13,120 (61.5%) | **≥ 100 (used)** | **11,945 (56.0%)** |
| **≥ 200 (used)** | **12,845 (60.2%)** | ≥ 200 | 9,729 (45.6%) |
| ≥ 500 | 11,250 (52.8%) | ≥ 500 | 4,293 (20.1%) |
| ≥ 1000 | 8,152 (38.2%) | ≥ 1000 | 960 (4.5%) |

At the chosen thresholds, ~13,200 windows contain any coding sequence at all and 12,845 (0D) /
11,945 (4D) survive, so the filter removes only 2.5% / 9.5% of gene-containing windows rather
than reshaping the sample. Going higher trades windows for precision quickly, especially for 4D:
`--min-sites-4d 500` would discard four fifths of windows and strongly favour the most gene-dense
ones.

The result is not sensitive to the cutoff. Spearman ρ for maize π<sub>4D</sub> is 0.336
(≥ 20 sites), 0.349 (≥ 100), 0.367 (≥ 200) and 0.408 (≥ 500); the modest increase is the
expected attenuation of noise, not a change in sign or story. The same holds for 0D
(0.266 → 0.325 across the same span) and for the 0D/4D ratio (−0.015 → −0.024 in maize).

## Results

12,845 windows pass for 0D, 11,945 for 4D, and 11,078–11,793 per population are paired for the
ratio (both classes pass and π<sub>4D</sub> > 0).

**π at 0-fold degenerate sites** (nonsynonymous proxy):

| Population | π<sub>0D</sub> | Spearman ρ | 95% CI | Partial ρ |
|------------|----------------|------------|--------|-----------|
| maize | 0.0034 | 0.271 | 0.238–0.306 | 0.187 |
| *mexicana* | 0.0053 | 0.205 | 0.171–0.235 | 0.147 |
| *parviglumis* | 0.0054 | 0.187 | 0.154–0.219 | 0.140 |

**π at 4-fold degenerate sites** (synonymous proxy):

| Population | π<sub>4D</sub> | Spearman ρ | 95% CI | Partial ρ |
|------------|----------------|------------|--------|-----------|
| maize | 0.0108 | 0.349 | 0.309–0.383 | 0.270 |
| *mexicana* | 0.0146 | 0.301 | 0.258–0.340 | 0.234 |
| *parviglumis* | 0.0158 | 0.293 | 0.254–0.330 | 0.243 |

π<sub>4D</sub> is close to the all-sites estimate for each population (0.0124 / 0.0205 / 0.0226
genome-wide), while π<sub>0D</sub> is roughly a third of it — the expected footprint of purifying
selection on replacement sites. The domestication bottleneck is visible in both classes: maize
carries about two thirds of teosinte diversity at 4D sites.

π by recombination decile, lowest (median 0.05 cM/Mb) to highest (median 3.31 cM/Mb):

| Population | 0D lowest | 0D highest | Fold | 4D lowest | 4D highest | Fold |
|------------|-----------|------------|------|-----------|------------|------|
| maize | 0.00263 | 0.00411 | 1.56× | 0.00636 | 0.01425 | 2.24× |
| *mexicana* | 0.00327 | 0.00649 | 1.98× | 0.00876 | 0.01862 | 2.13× |
| *parviglumis* | 0.00379 | 0.00605 | 1.60× | 0.01004 | 0.01932 | 1.92× |

![π at 0D and 4D sites vs recombination rate](results/pi_coding_vs_recombination.png)

Both classes reproduce the all-sites shape — a steep rise across the low-recombination deciles
that saturates above roughly 0.5 cM/Mb. The correlation is consistently *stronger at 4D than at
0D* in every population, which is what linked selection predicts: neutral diversity is free to
track the local effective population size, while diversity at constrained sites is held down by
direct purifying selection regardless of the local rate.

**The 0D/4D ratio.** Genome-wide (over paired windows) π<sub>0D</sub>/π<sub>4D</sub> is 0.310 in
maize, 0.358 in *mexicana* and 0.338 in *parviglumis* — maize lowest, though with only 8
haplotypes this difference should not be pushed hard.

| Population | π<sub>0D</sub>/π<sub>4D</sub> | Spearman ρ vs rate | 95% CI | Partial ρ | Lowest decile | Highest decile |
|------------|-------------------------------|--------------------|--------|-----------|---------------|----------------|
| maize | 0.310 | −0.021 | −0.041 to −0.001 | −0.038 | 0.388 | 0.285 |
| *mexicana* | 0.358 | −0.031 | −0.056 to −0.005 | −0.033 | 0.355 | 0.354 |
| *parviglumis* | 0.338 | −0.050 | −0.074 to −0.026 | −0.055 | 0.362 | 0.315 |

![π0D/π4D vs recombination rate](results/pi_coding_vs_recombination_ratio.png)

The ratio does decline with recombination, in the predicted direction, in all three populations,
and the bootstrap interval excludes zero in all three. But the effect is **an order of magnitude
weaker than the raw π–rate correlations** (|ρ| ≈ 0.02–0.05 against 0.19–0.35) and the decile
trend is non-monotonic, with a conspicuous high point in the lowest deciles and a dip near
1.2 cM/Mb. Read honestly: this is a weak signal that is consistent with more efficient purifying
selection in high-recombination regions, not a demonstration of it. Most of what drives π up with
recombination is shared between 0D and 4D sites and therefore cancels in the ratio, which is
precisely why the ratio is the interpretable quantity and also why it is small.

## Caveats

**Autocorrelation.** As in the all-sites analysis, adjacent windows are not independent; all
intervals come from a moving-block bootstrap over contiguous 5 Mb blocks within chromosomes, and
no asymptotic p-value is reported.

**Gene density is the central confounder here.** 0D and 4D sites exist only in genes, gene
density rises toward the chromosome arms, and so does recombination rate. Controlling for gene
density and centromere distance shrinks the correlations by roughly a third (maize 4D
0.349 → 0.270; maize 0D 0.271 → 0.187) — less than the halving seen in the all-sites analysis,
but the residual association is still not free of shared chromosome-scale structure. For the
ratio, conditioning barely changes anything (maize −0.021 → −0.038), as expected for a quantity
in which gene-density effects largely cancel.

**Ascertainment runs the other way from the all-sites analysis.** There, callable-site density
was *negatively* correlated with rate (ρ ≈ −0.28). Here it is *positively* correlated
(ρ ≈ +0.36 for both classes), because coding sequence is both more alignable and more abundant
in high-recombination regions. That is not conservative: any tendency for windows with more
callable coding sites to yield better-estimated (or systematically different) π is aligned with
the rate axis rather than opposed to it. The ratio is again the more robust statistic, since 0D
and 4D callability track each other almost exactly.

**Reference-codon degeneracy.** Classes are assigned from the B73 reference codon only
(see [above](#why-0d4d-rather-than-literal-synonymousnonsynonymous)), and B73 is a maize inbred,
so the annotation is very slightly closer to the maize samples than to the teosintes. Any
resulting bias affects all three populations in the same direction and is small relative to the
effects reported.

**Annotation quality.** Everything rests on the Zm00001eb.1 gene models. Misannotated genes,
pseudogenes annotated as protein-coding, and transposon-derived ORFs all inject sites whose
degeneracy labels are meaningless; no expression or conservation filter was applied.

**Sample size.** 8, 12 and 9 haplotypes over a few hundred sites per window make per-window π —
and especially the per-window ratio — very noisy. The scatter panel of the ratio figure shows
that spread plainly. Only the aggregate trends should be interpreted; individual windows should
not.
