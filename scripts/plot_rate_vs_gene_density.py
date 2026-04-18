#!/usr/bin/env python3
"""Plot recombination rate (cM/Mb) vs gene density in 100 kb windows."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent
WIN = 100_000
CHROMS = [f"chr{i}" for i in range(1, 11)]


def load_finemap():
    fm = pd.read_csv(
        ROOT / "data/finemap_v5.bed", sep="\t", header=None,
        names=["chr", "start", "end", "cM_start", "cM_end", "cM_per_Mb"],
    )
    fm["chr"] = fm["chr"].str.lower()
    return fm[fm["chr"].isin(CHROMS)].sort_values(["chr", "start"])


def load_genes():
    gff = pd.read_csv(
        ROOT / "data/v5.genes.gff3", sep="\t", header=None, comment="#",
        usecols=[0, 2, 3, 4], names=["chr", "feature", "start", "end"],
    )
    gff = gff[gff["feature"] == "gene"].copy()
    gff["chr"] = gff["chr"].str.lower()
    gff = gff[gff["chr"].isin(CHROMS)]
    gff["mid"] = (gff["start"] + gff["end"]) // 2
    return gff


def load_chr_lengths():
    fai = pd.read_csv(
        ROOT / "data/v5.fa.gz.fai", sep="\t", header=None,
        usecols=[0, 1], names=["chr", "length"],
    )
    fai["chr"] = fai["chr"].str.lower()
    return fai[fai["chr"].isin(CHROMS)].set_index("chr")["length"]


def compute_windows(fm, genes, chr_lengths):
    records = []
    for chrom in CHROMS:
        if chrom not in chr_lengths.index:
            continue
        chrlen = chr_lengths[chrom]
        fm_c = fm[fm["chr"] == chrom]
        gene_mids = genes.loc[genes["chr"] == chrom, "mid"].values

        for win_start in range(0, int(chrlen), WIN):
            win_end = min(win_start + WIN, int(chrlen))

            segs = fm_c[(fm_c["end"] > win_start) & (fm_c["start"] < win_end)]
            if segs.empty:
                continue
            overlap = (
                np.minimum(segs["end"].values, win_end)
                - np.maximum(segs["start"].values, win_start)
            ).clip(min=0)
            total_overlap = overlap.sum()
            if total_overlap == 0:
                continue
            avg_rate = float((segs["cM_per_Mb"].values * overlap).sum() / total_overlap)

            n_genes = int(np.sum((gene_mids >= win_start) & (gene_mids < win_end)))
            gene_density = n_genes / (WIN / 1e6)

            records.append({
                "chr": chrom,
                "win_start": win_start,
                "cM_per_Mb": avg_rate,
                "gene_density": gene_density,
            })
    return pd.DataFrame(records)


def binned_medians(x, y, n_bins=20):
    bins = np.percentile(x, np.linspace(0, 100, n_bins + 1))
    bins = np.unique(bins)
    idx = np.digitize(x, bins) - 1
    med_x, med_y = [], []
    for b in range(len(bins) - 1):
        mask = idx == b
        if mask.sum() > 0:
            med_x.append(np.median(x[mask]))
            med_y.append(np.median(y[mask]))
    return np.array(med_x), np.array(med_y)


def plot(df):
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(
        df["cM_per_Mb"], df["gene_density"],
        s=4, alpha=0.3, linewidths=0, color="#2166ac", zorder=2,
    )
    med_x, med_y = binned_medians(
        df["cM_per_Mb"].values, df["gene_density"].values
    )
    ax.plot(med_x, med_y, color="#d62728", linewidth=2,
            label="Binned median", zorder=3)
    ax.set_xlabel("Recombination rate (cM/Mb)", fontsize=11)
    ax.set_ylabel("Gene density (genes / Mb)", fontsize=11)
    ax.set_title(
        "Gene density vs recombination rate\n(100 kb windows, chr1–10)", fontsize=11
    )
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = ROOT / "results/rate_vs_gene_density.png"
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")


def main():
    print("Loading data...")
    fm = load_finemap()
    genes = load_genes()
    chr_lengths = load_chr_lengths()
    print("Computing windows...")
    df = compute_windows(fm, genes, chr_lengths)
    print(f"  {len(df)} windows")
    print(df[["cM_per_Mb", "gene_density"]].describe().round(3))
    plot(df)


if __name__ == "__main__":
    main()
