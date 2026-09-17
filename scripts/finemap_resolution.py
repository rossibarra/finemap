#!/usr/bin/env python3
"""Diagnose the fine-scale resolution limit of the FineMap interval-density map.

build_finemap.py assigns each crossover a weight of 1/(end - start) and spreads it
uniformly across its interval.  Per-bp weight density therefore scales as 1/width^2,
so a handful of unusually narrow crossover intervals dominate all fine-scale structure
in the map while the typical interval contributes a broad, flat smear.

This script quantifies that: the width distribution of the source crossover intervals,
the 1/width^2 weighting, and -- if a hotspot BED is supplied -- the share of each
apparent hotspot's crossover weight that comes from atypically narrow intervals.
"""

import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NARROW = 10_000  # "narrow" crossover interval threshold, bp


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jri", required=True, help="crossover intervals BED (data/jri_v5.bed)")
    ap.add_argument("--hotspots", help="hotspot BED to diagnose")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    jri = pd.read_csv(args.jri, sep="\t", header=None,
                      names=["chrom", "start", "end", "id", "src"])
    jri["chrom"] = jri.chrom.str.replace("Chr", "chr", regex=False)
    jri = jri[jri.end > jri.start].copy()
    jri["w"] = jri.end - jri.start
    jri["wt"] = 1.0 / jri.w

    med = float(jri.w.median())
    print(f"crossover intervals        {len(jri):,}")
    print(f"  median width             {med:,.0f} bp")
    print(f"  mean width               {jri.w.mean():,.0f} bp")
    for q in (5, 25, 50, 75, 95, 99):
        print(f"  p{q:<2d}                      {np.percentile(jri.w, q):,.0f} bp")
    print(f"  fraction < {NARROW:,} bp      {(jri.w < NARROW).mean():.1%}")
    print(f"  fraction > 100 kb        {(jri.w > 100_000).mean():.1%}")

    frac = None
    if args.hotspots:
        hs = pd.read_csv(args.hotspots, sep="\t", header=None,
                         usecols=[0, 1, 2], names=["chrom", "start", "end"])
        rows = []
        for chrom, g in hs.groupby("chrom"):
            j = jri[jri.chrom == chrom]
            js, je = j.start.to_numpy(), j.end.to_numpy()
            wt, ww = j.wt.to_numpy(), j.w.to_numpy()
            for s, e in zip(g.start, g.end):
                m = (js < e) & (je > s)
                if not m.any():
                    continue
                # weight actually deposited inside the hotspot by each source interval
                overlap = np.minimum(je[m], e) - np.maximum(js[m], s)
                contrib = wt[m] * overlap
                if contrib.sum() > 0:
                    rows.append(contrib[ww[m] < NARROW].sum() / contrib.sum())
        frac = np.array(rows)
        print(f"\nhotspots examined          {len(frac):,}")
        print(f"  share of crossover weight from intervals < {NARROW:,} bp:")
        print(f"    median                 {np.median(frac):.1%}")
        print(f"    mean                   {frac.mean():.1%}")
        print(f"    hotspots above 50%     {(frac > 0.5).mean():.0%}")

    n_panels = 3 if frac is not None else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(4.3 * n_panels, 3.9))

    ax = axes[0]
    ax.hist(np.log10(jri.w), bins=70, color="#3B6FA0", alpha=0.85)
    ax.axvline(np.log10(med), color="#B8860B", lw=1.6,
               label=f"median {med/1000:.0f} kb")
    ax.axvline(np.log10(NARROW), color="#555555", lw=1.2, ls="--",
               label=f"{NARROW//1000} kb")
    ax.set_xlabel("crossover interval width (log₁₀ bp)")
    ax.set_ylabel("intervals")
    ax.set_title("Source crossover intervals\nare wide")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    w = np.logspace(2, 7, 200)
    ax.loglog(w, 1 / w**2, color="#2E7D5B", lw=1.8)
    ax.axvline(med, color="#B8860B", lw=1.6)
    ax.annotate(f"a 100 bp interval carries\n{(1/100**2)/(1/med**2):,.0f}× the per-bp weight\nof a median interval",
                xy=(100, 1e-4), xytext=(2.5e3, 2e-5), fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.9, color="#555555"))
    ax.set_xlabel("interval width (bp)")
    ax.set_ylabel("weight per bp  (1/width²)")
    ax.set_title("Narrow intervals dominate\nfine-scale structure")

    if frac is not None:
        ax = axes[2]
        ax.hist(100 * frac, bins=30, color="#B8860B", alpha=0.85)
        ax.axvline(100 * np.median(frac), color="#333333", lw=1.6,
                   label=f"median {np.median(frac):.0%}")
        ax.set_xlabel(f"% of hotspot CO weight from intervals < {NARROW//1000} kb")
        ax.set_ylabel("hotspots")
        ax.set_title("Apparent hotspots are built\nfrom the narrow tail")
        ax.legend(frameon=False, fontsize=8)

    for a in axes:
        style(a)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
