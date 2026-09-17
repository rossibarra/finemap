#!/usr/bin/env python3
"""Compare pi at all sites, 4-fold and 0-fold degenerate sites against recombination rate.

The three site classes each need their own callable-site threshold, so their separate
analyses rest on different window sets.  This script restricts all three to the windows
that pass every filter, making the comparison between classes exact: same windows, same
recombination range, differing only in which sites are counted.

Top row: pi by recombination decile, as sum(count_diffs)/sum(count_comparisons) within
each decile, which keeps zero-pi windows in the estimate where they belong.
Bottom row: every window, log-log, with an ordinary least-squares fit of log10(pi) on
log10(rate).  Windows with pi = 0 cannot be shown on a log axis and are dropped from the
scatter and the fit only; the fraction dropped is printed and annotated.
"""

import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

POP_COLORS = {"maize": "#B8860B", "parviglumis": "#2E7D5B", "mexicana": "#3B6FA0"}
POP_ORDER = ["maize", "mexicana", "parviglumis"]
CLASSES = [("all sites", "pi_finemap_windows.tsv", 5000),
           ("4D", "pi_4D_finemap_windows.tsv", 100),
           ("0D", "pi_0D_finemap_windows.tsv", 200)]
TITLES = {"all sites": "All sites", "4D": "4-fold degenerate (synonymous)",
          "0D": "0-fold degenerate (nonsynonymous)"}


def spearman(x, y):
    rx = np.array(pd.Series(x).rank(), dtype=float)
    ry = np.array(pd.Series(y).rank(), dtype=float)
    rx -= rx.mean(); ry -= ry.mean()
    return float((rx * ry).sum() / np.sqrt((rx * rx).sum() * (ry * ry).sum()))


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pi-dir", default="data/pixy")
    ap.add_argument("--bed", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    bed = pd.read_csv(args.bed, sep="\t", header=None,
                      names=["chrom", "start", "end", "cs", "ce", "rate"])
    bed["chrom"] = bed.chrom.str.replace("Chr", "chr", regex=False)

    tabs, keysets = {}, []
    for cls, fname, min_sites in CLASSES:
        d = pd.read_csv(f"{args.pi_dir}/{fname}", sep="\t").rename(columns={"chromosome": "chrom"})
        d["start"] = d.window_pos_1 - 1
        d = d[d.avg_pi.notna() & (d.no_sites >= min_sites)]
        tabs[cls] = d
        sub = d[d["pop"] == POP_ORDER[0]][["chrom", "start"]]
        keysets.append(set(map(tuple, sub.to_numpy())))

    common = pd.DataFrame(sorted(set.intersection(*keysets)), columns=["chrom", "start"])
    print(f"windows passing all three filters: {len(common):,}\n")

    data, rows = {}, []
    for cls in tabs:
        m = (tabs[cls].merge(common, on=["chrom", "start"])
                      .merge(bed[["chrom", "start", "rate"]], on=["chrom", "start"]))
        data[cls] = m
        for pop, g in m.groupby("pop"):
            pos = g[g.avg_pi > 0]
            slope, intercept = np.polyfit(np.log10(pos.rate), np.log10(pos.avg_pi), 1)
            rows.append({"class": cls, "pop": pop, "n_windows": len(g),
                         "pi": g.count_diffs.sum() / g.count_comparisons.sum(),
                         "spearman_rho": spearman(g.rate, g.avg_pi),
                         "loglog_slope": slope,
                         "pct_pi_zero": 100 * (g.avg_pi <= 0).mean()})

    summary = pd.DataFrame(rows)
    summary["class"] = pd.Categorical(summary["class"], [c[0] for c in CLASSES], ordered=True)
    summary["pop"] = pd.Categorical(summary["pop"], POP_ORDER, ordered=True)
    summary = summary.sort_values(["class", "pop"])
    summary.to_csv(f"{args.out_prefix}_summary.tsv", sep="\t", index=False, float_format="%.4f")
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.6))
    for j, (cls, _, _) in enumerate(CLASSES):
        m = data[cls]

        ax = axes[0, j]
        for pop in POP_ORDER:
            g = m[m["pop"] == pop]
            q = pd.qcut(g.rate.rank(method="first"), 10, labels=False)
            b = g.groupby(q).apply(lambda x: pd.Series({
                "rate": x.rate.median(),
                "pi": x.count_diffs.sum() / x.count_comparisons.sum(),
                "se": x.avg_pi.std() / np.sqrt(len(x))}), include_groups=False)
            ax.errorbar(b.rate, b.pi, yerr=b.se, marker="o", ms=4, lw=1.5, capsize=2,
                        color=POP_COLORS[pop], label=pop)
        ax.set_xscale("log")
        ax.set_xlabel("recombination rate (cM/Mb), decile median")
        ax.set_ylabel(r"$\pi$" if j == 0 else "")
        ax.set_title(TITLES[cls], fontsize=11)
        if j == 0:
            ax.legend(frameon=False, fontsize=9)

        ax = axes[1, j]
        for pop in POP_ORDER:
            g = m[(m["pop"] == pop) & (m.avg_pi > 0)]
            ax.scatter(g.rate, g.avg_pi, s=2, alpha=0.10, color=POP_COLORS[pop],
                       rasterized=True, linewidths=0)
        xs = np.logspace(np.log10(m.rate.min()), np.log10(m.rate.max()), 100)
        for pop in POP_ORDER:
            g = m[(m["pop"] == pop) & (m.avg_pi > 0)]
            s, b0 = np.polyfit(np.log10(g.rate), np.log10(g.avg_pi), 1)
            ax.plot(xs, 10 ** (b0 + s * np.log10(xs)), color=POP_COLORS[pop], lw=2,
                    label=f"{pop}  slope {s:.2f}")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("recombination rate (cM/Mb)")
        ax.set_ylabel(r"$\pi$" if j == 0 else "")
        dropped = 100 * (m.avg_pi <= 0).mean()
        if dropped > 0:
            ax.text(0.03, 0.04, f"{dropped:.1f}% of windows have π = 0\nand are omitted",
                    transform=ax.transAxes, fontsize=7.5, color="#666666", va="bottom")
        ax.legend(frameon=False, fontsize=8, loc="upper left")

    for ax in axes.ravel():
        style(ax)
    fig.suptitle(f"Diversity vs recombination rate by site class "
                 f"({len(common):,} matched 100 kb windows)", y=0.995, fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.png", dpi=200)
    print(f"\nwrote {args.out_prefix}.png, {args.out_prefix}_summary.tsv")


if __name__ == "__main__":
    main()
