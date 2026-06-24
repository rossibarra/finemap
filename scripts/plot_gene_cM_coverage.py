#!/usr/bin/env python3
"""Per-chromosome bar plot: % physical bp vs % cM within gene flanks."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import argparse

ROOT = Path(__file__).parent.parent
CHROMS = [f"chr{i}" for i in range(1, 11)]


def load_chr_lengths():
    fai = pd.read_csv(
        ROOT / "data/v5.fa.gz.fai", sep="\t", header=None,
        usecols=[0, 1], names=["chr", "length"],
    )
    fai["chr"] = fai["chr"].str.lower()
    return fai[fai["chr"].isin(CHROMS)].set_index("chr")["length"]


def load_genes():
    gff = pd.read_csv(
        ROOT / "data/v5.genes.gff3", sep="\t", header=None, comment="#",
        usecols=[0, 2, 3, 4], names=["chr", "feature", "start", "end"],
    )
    gff = gff[gff["feature"] == "gene"].copy()
    gff["chr"] = gff["chr"].str.lower()
    return gff[gff["chr"].isin(CHROMS)]


def merge_intervals(starts, ends):
    """Merge sorted overlapping intervals; returns array of (start, end)."""
    if len(starts) == 0:
        return np.empty((0, 2), dtype=np.int64)
    order = np.argsort(starts)
    starts, ends = starts[order], ends[order]
    merged = [[starts[0], ends[0]]]
    for s, e in zip(starts[1:], ends[1:]):
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return np.array(merged, dtype=np.int64)


def covered_bp(merged, chr_len):
    """Total bp covered by merged intervals, clamped to [0, chr_len]."""
    clipped = np.clip(merged, 0, chr_len)
    return int(np.sum(np.maximum(0, clipped[:, 1] - clipped[:, 0])))


def overlap_with_merged(seg_starts, seg_ends, merged):
    """
    For each segment [seg_starts[i], seg_ends[i]], compute bp overlapping
    with any interval in `merged`. Returns array of overlap lengths.
    """
    overlaps = np.zeros(len(seg_starts), dtype=np.float64)
    m_starts = merged[:, 0]
    m_ends = merged[:, 1]
    for i, (ss, se) in enumerate(zip(seg_starts, seg_ends)):
        # find merged intervals that could overlap [ss, se]
        lo = np.searchsorted(m_ends, ss, side="right")
        hi = np.searchsorted(m_starts, se, side="left")
        if hi <= lo:
            continue
        ov = np.sum(
            np.minimum(m_ends[lo:hi], se) - np.maximum(m_starts[lo:hi], ss)
        )
        overlaps[i] = max(0, ov)
    return overlaps


def flank_label(flank):
    return f"{flank // 1000} kb" if flank >= 1000 and flank % 1000 == 0 else f"{flank} bp"


def compute_coverage(flank, chr_lengths, genes, finemap):
    pct_bp, pct_cM = [], []

    for chrom in CHROMS:
        chr_len = int(chr_lengths[chrom])
        g = genes[genes["chr"] == chrom]

        starts = np.maximum(0, g["start"].values - flank).astype(np.int64)
        ends = np.minimum(chr_len, g["end"].values + flank).astype(np.int64)
        merged = merge_intervals(starts, ends)

        # % physical
        bp_cov = covered_bp(merged, chr_len)
        pct_bp.append(100.0 * bp_cov / chr_len)

        # % cM
        fm = finemap[finemap["chr"] == chrom].sort_values("start")
        seg_s = fm["start"].values.astype(np.int64)
        seg_e = fm["end"].values.astype(np.int64)
        seg_cM = (fm["cM_end"] - fm["cM_start"]).values
        total_cM = seg_cM.sum()

        if total_cM == 0 or len(merged) == 0:
            pct_cM.append(0.0)
            continue

        seg_len = (seg_e - seg_s).astype(np.float64)
        ov = overlap_with_merged(seg_s, seg_e, merged)
        cM_in = np.sum(seg_cM * np.where(seg_len > 0, ov / seg_len, 0.0))
        pct_cM.append(100.0 * cM_in / total_cM)

    return pct_bp, pct_cM


def plot(results):
    x = np.arange(len(CHROMS))
    n_flanks = len(results)
    fig, axes = plt.subplots(
        n_flanks, 1, figsize=(10, 4.2 * n_flanks), sharex=True, squeeze=False
    )
    axes = axes.flatten()

    for ax, (flank, pct_bp, pct_cM) in zip(axes, results):
        w = 0.35
        label = flank_label(flank)
        ax.bar(x - w / 2, pct_bp, width=w, label="% physical bp", color="#4393c3")
        ax.bar(x + w / 2, pct_cM, width=w, label="% cM", color="#d6604d")
        ax.set_ylabel(f"% within gene ± {label}")
        ax.set_title(f"Gene ± {label} coverage: physical vs genetic (cM)", fontsize=11)
        ax.legend(fontsize=9)

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels([f"Chr{i}" for i in range(1, 11)])
    fig.tight_layout()
    out = ROOT / "results/gene_cM_coverage.png"
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")


def print_table(results):
    headers = ["Chr"]
    for flank, _, _ in results:
        label = flank_label(flank)
        headers.append(f"%bp (±{label})")
    for flank, _, _ in results:
        label = flank_label(flank)
        headers.append(f"%cM (±{label})")
    print("\n" + "  ".join(f"{h:>12}" for h in headers))
    for idx, chrom in enumerate(CHROMS):
        row = [chrom]
        for _, pct_bp, pct_cM in results:
            row.append(f"{pct_bp[idx]:.1f}")
        for _, pct_bp, pct_cM in results:
            row.append(f"{pct_cM[idx]:.1f}")
        print("  ".join(f"{value:>12}" for value in row))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--flank",
        type=int,
        nargs="+",
        default=[1_000],
        help="One or more flank sizes in bp around each gene (default: 1000)",
    )
    args = parser.parse_args()

    chr_lengths = load_chr_lengths()
    genes = load_genes()
    finemap = pd.read_csv(
        ROOT / "data/finemap_v5.bed", sep="\t", header=None,
        names=["chr", "start", "end", "cM_start", "cM_end", "cM_per_Mb"],
    )
    finemap["chr"] = finemap["chr"].str.lower()
    finemap = finemap[finemap["chr"].isin(CHROMS)]

    results = []
    for flank in args.flank:
        pct_bp, pct_cM = compute_coverage(flank, chr_lengths, genes, finemap)
        results.append((flank, pct_bp, pct_cM))
    plot(results)
    print_table(results)


if __name__ == "__main__":
    main()
