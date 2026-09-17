#!/usr/bin/env python3
"""Nucleotide diversity as a function of distance to the nearest recombination hotspot.

Takes a pixy-format pi table on fine (2 kb) windows, a hotspot BED, the smoothed
hierarchical map (local background rate).  For every window it computes the distance
from the window midpoint to the nearest hotspot on the same chromosome, bins windows on a log distance scale, and reports pi per bin per population
as sum(count_diffs) / sum(count_comparisons) -- never the mean of per-window pi, which
would let sparse, noisy windows dominate.

Hotspots sit in high-recombination distal sequence, which raises pi on its own.  The
script therefore also reports background rate as a function of distance, repeats the
distance binning inside quartiles of local rate, and reports a partial Spearman of pi
against log distance controlling for local rate.
"""

import argparse
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", message="All-NaN slice encountered")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

POP_COLORS = {"maize": "#B8860B", "parviglumis": "#2E7D5B", "mexicana": "#3B6FA0"}
POP_ORDER = ["maize", "mexicana", "parviglumis"]
BED_COLUMNS = ["chrom", "start", "end", "cm_start", "cm_end", "rate"]

# distance bin edges in bp; the first bin is "overlaps a hotspot" (distance 0)
EDGES = [0, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, np.inf]
LABELS = ["0", "0-1k", "1-2k", "2-5k", "5-10k", "10-20k", "20-50k",
          "50-100k", "100-200k", "200-500k", ">500k"]


def spearman(x, y):
    """Spearman rho without scipy: Pearson on ranks."""
    rx = np.array(pd.Series(x).rank(), dtype=float)
    ry = np.array(pd.Series(y).rank(), dtype=float)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / denom) if denom > 0 else np.nan


def partial_spearman(x, y, covariates):
    """Spearman of x and y after regressing rank-transformed covariates out of both."""
    rx = np.array(pd.Series(x).rank(), dtype=float)
    ry = np.array(pd.Series(y).rank(), dtype=float)
    cov = np.column_stack([np.array(pd.Series(c).rank(), dtype=float) for c in covariates])
    design = np.column_stack([np.ones(len(rx)), cov])
    rx_r = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    ry_r = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    return spearman(rx_r, ry_r)


def hotspot_distance(chrom, start, end, hot):
    """Distance from each window midpoint to the nearest hotspot on the same chromosome.

    Zero when the window overlaps a hotspot.  Uses searchsorted per chromosome.
    """
    mid = (start + end) / 2.0
    dist = np.full(len(mid), np.inf)
    for c, h in hot.groupby("chrom", sort=False):
        m = chrom == c
        if not m.any():
            continue
        hs = np.sort(h.start.to_numpy())
        order = np.argsort(h.start.to_numpy())
        he = h.end.to_numpy()[order]
        mm = mid[m]
        i = np.searchsorted(hs, mm, side="right")      # hotspots strictly left: i-1
        left = np.where(i > 0, mm - he[np.clip(i - 1, 0, len(he) - 1)], np.inf)
        left = np.where(i > 0, np.maximum(left, 0.0), np.inf)
        right = np.where(i < len(hs), hs[np.clip(i, 0, len(hs) - 1)] - mm, np.inf)
        right = np.where(i < len(hs), np.maximum(right, 0.0), np.inf)
        d = np.minimum(left, right)
        # a window whose midpoint falls inside a hotspot gets 0 from the above only if
        # the hotspot starts at/below mid and ends above it; force that case explicitly
        inside = (i > 0) & (he[np.clip(i - 1, 0, len(he) - 1)] > mm)
        d = np.where(inside, 0.0, d)
        # also zero when the window interval itself overlaps a hotspot
        s, e = start[m], end[m]
        j = np.searchsorted(hs, e, side="left")
        ov = (j > 0) & (he[np.clip(j - 1, 0, len(he) - 1)] > s)
        d = np.where(ov, 0.0, d)
        dist[m] = d
    return dist


def local_rate(chrom, mid, hier):
    """Rate from the uniform hierarchical map at each window midpoint."""
    size = int((hier.end - hier.start).mode().iloc[0])
    lookup = {}
    for c, g in hier.groupby("chrom", sort=False):
        idx = (g.start.to_numpy() // size).astype(int)
        arr = np.full(idx.max() + 1, np.nan)
        arr[idx] = g.rate.to_numpy()
        lookup[c] = arr
    out = np.full(len(mid), np.nan)
    for c, arr in lookup.items():
        m = chrom == c
        if not m.any():
            continue
        k = np.clip((mid[m] // size).astype(int), 0, len(arr) - 1)
        out[m] = arr[k]
    return out


def ratio_pi(diffs, pairs):
    return float(diffs.sum() / pairs.sum()) if pairs.sum() > 0 else np.nan


def block_bootstrap_bins(df, block_id, n_bins, n_boot, seed):
    """Percentile CI for ratio-of-sums pi in each distance bin.

    Resamples contiguous 1 Mb blocks of windows with replacement, which preserves the
    along-chromosome autocorrelation that makes naive per-window SEs too small.
    """
    rng = np.random.default_rng(seed)
    blocks, binv = pd.factorize(block_id)
    nb = blocks.max() + 1
    D = np.zeros((nb, n_bins))
    P = np.zeros((nb, n_bins))
    np.add.at(D, (blocks, df.bin_idx.to_numpy()), df.count_diffs.to_numpy())
    np.add.at(P, (blocks, df.bin_idx.to_numpy()), df.count_comparisons.to_numpy())
    reps = np.empty((n_boot, n_bins))
    for b in range(n_boot):
        pick = rng.integers(0, nb, size=nb)
        d = D[pick].sum(axis=0)
        p = P[pick].sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            reps[b] = np.where(p > 0, d / p, np.nan)
    lo = np.nanquantile(reps, 0.025, axis=0)
    hi = np.nanquantile(reps, 0.975, axis=0)
    return lo, hi


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pi", required=True, help="pixy-format pi table on fine windows")
    ap.add_argument("--hotspots", required=True, help="hotspot BED (chrom start end ...)")
    ap.add_argument("--hierarchical", required=True, help="uniform hierarchical map BED")
    ap.add_argument("--min-sites", type=int, default=200,
                    help="drop windows with fewer callable sites (default 200)")
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--block-mb", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    hot = pd.read_csv(args.hotspots, sep="\t", header=None,
                      names=["chrom", "start", "end", "cm", "rate", "n_intervals"])
    hot["chrom"] = hot.chrom.str.replace("Chr", "chr", regex=False)
    hier = pd.read_csv(args.hierarchical, sep="\t", header=None, names=BED_COLUMNS)
    hier["chrom"] = hier.chrom.str.replace("Chr", "chr", regex=False)

    pi = pd.read_csv(args.pi, sep="\t")
    pi = pi.rename(columns={"chromosome": "chrom"})
    pi["start"] = pi.window_pos_1 - 1
    pi["end"] = pi.window_pos_2

    n_win_total = len(pi) // pi["pop"].nunique()
    print(f"fine windows: {n_win_total:,} per population")
    med = pi.groupby("pop").no_sites.median()
    print(f"median callable sites per window: {med.iloc[0]:,.0f}")
    print("callable-site percentiles:",
          {f"p{q}": int(np.percentile(pi[pi["pop"] == POP_ORDER[0]].no_sites, q))
           for q in (10, 25, 50, 75, 90)})

    keep = pi.avg_pi.notna() & (pi.no_sites >= args.min_sites)
    df = pi[keep].copy()
    print(f"windows passing >= {args.min_sites} callable sites: "
          f"{len(df) // df['pop'].nunique():,} of {n_win_total:,} "
          f"({100 * len(df) / len(pi):.1f}%)")

    chrom = df.chrom.to_numpy()
    mid = ((df.start.to_numpy() + df.end.to_numpy()) / 2.0)
    df["dist"] = hotspot_distance(chrom, df.start.to_numpy(), df.end.to_numpy(), hot)
    df["local_rate"] = local_rate(chrom, mid, hier)
    df = df[np.isfinite(df.dist) & df.local_rate.notna()].copy()

    # bin 0 is "window overlaps a hotspot"; bins 1.. are the log-spaced distance classes
    hi_edges = np.array(EDGES[1:])
    idx = np.searchsorted(hi_edges, df.dist.to_numpy(), side="left") + 1
    idx[df.dist.to_numpy() == 0] = 0
    df["bin_idx"] = idx
    n_bins = len(LABELS)
    block = df.chrom + ":" + (df.start // int(args.block_mb * 1e6)).astype(str)

    # ---------------- main table: pi by distance bin ----------------
    rows = []
    for pop in POP_ORDER:
        g = df[df["pop"] == pop]
        lo, hi = block_bootstrap_bins(g, block[g.index], n_bins, args.n_boot, args.seed)
        for b in range(n_bins):
            sub = g[g.bin_idx == b]
            if len(sub) == 0:
                continue
            rows.append({
                "pop": pop, "bin": LABELS[b], "bin_idx": b, "n_windows": len(sub),
                "median_dist_bp": float(sub.dist.median()),
                "pi": ratio_pi(sub.count_diffs, sub.count_comparisons),
                "ci_low": lo[b], "ci_high": hi[b],
                "median_local_rate": float(sub.local_rate.median()),
                "median_callable": float(sub.no_sites.median()),
            })
    binned = pd.DataFrame(rows)
    binned.to_csv(f"{args.out_prefix}_bins.tsv", sep="\t", index=False, float_format="%.6g")

    print("\n=== pi by distance to nearest hotspot ===")
    for pop in POP_ORDER:
        t = binned[binned["pop"] == pop]
        print(f"\n{pop}")
        print(t[["bin", "n_windows", "pi", "ci_low", "ci_high",
                 "median_local_rate", "median_callable"]]
              .to_string(index=False, float_format=lambda v: f"{v:.5g}"))

    # ---------------- stratified by local recombination quartile ----------------
    ref = df[df["pop"] == POP_ORDER[0]]
    qedges = np.quantile(ref.local_rate, [0.25, 0.5, 0.75])
    df["rate_q"] = np.searchsorted(qedges, df.local_rate.to_numpy(), side="right")
    print("\nlocal-rate quartile breaks (cM/Mb): " +
          ", ".join(f"{v:.3f}" for v in qedges))

    srows = []
    for q in range(4):
        dq = df[df.rate_q == q]
        for pop in POP_ORDER:
            g = dq[dq["pop"] == pop]
            if len(g) == 0:
                continue
            lo, hi = block_bootstrap_bins(g, block[g.index], n_bins,
                                          max(200, args.n_boot // 2), args.seed)
            for b in range(n_bins):
                sub = g[g.bin_idx == b]
                if len(sub) < 20:
                    continue
                srows.append({"rate_q": q, "pop": pop, "bin": LABELS[b], "bin_idx": b,
                              "n_windows": len(sub),
                              "pi": ratio_pi(sub.count_diffs, sub.count_comparisons),
                              "ci_low": lo[b], "ci_high": hi[b],
                              "median_local_rate": float(sub.local_rate.median())})
    strat = pd.DataFrame(srows)
    strat.to_csv(f"{args.out_prefix}_stratified.tsv", sep="\t", index=False,
                 float_format="%.6g")

    # ---------------- correlations ----------------
    crows = []
    logd = np.log10(df.dist.to_numpy() + 1000.0)
    df["logd"] = logd
    for pop in POP_ORDER:
        g = df[df["pop"] == pop]
        rho = spearman(g.logd, g.avg_pi)
        prho = partial_spearman(g.logd, g.avg_pi, [g.local_rate])
        # delete-one-chromosome jackknife SE
        jk_full = []
        for c in sorted(g.chrom.unique()):
            sub = g[g.chrom != c]
            jk_full.append(spearman(sub.logd, sub.avg_pi))
        jk_full = np.array(jk_full)
        n = len(jk_full)
        jk_se = np.sqrt((n - 1) / n * ((jk_full - jk_full.mean()) ** 2).sum())
        crows.append({"pop": pop, "n_windows": len(g),
                      "rho_pi_vs_logdist": rho, "jackknife_se": jk_se,
                      "partial_rho_given_rate": prho,
                      "rho_logdist_vs_local_rate": spearman(g.logd, g.local_rate),
                      "rho_callable_vs_logdist": spearman(g.logd, g.no_sites),
                      "genome_pi": ratio_pi(g.count_diffs, g.count_comparisons)})
    corr = pd.DataFrame(crows)
    corr.to_csv(f"{args.out_prefix}_summary.tsv", sep="\t", index=False,
                float_format="%.4f")
    print("\n=== correlations (per 2 kb window) ===")
    print(corr.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    # ---------------- plots ----------------
    def xpos(b):
        """Log-scale x position in kb; the overlap bin sits left of the axis break."""
        if b == 0:
            return 0.32
        lo_, hi_ = EDGES[b - 1], EDGES[b]
        if not np.isfinite(hi_):
            return 800.0
        return np.sqrt(max(lo_, 500.0) * hi_) / 1e3

    xs = np.array([xpos(b) for b in range(n_bins)])

    def style(ax):
        ax.set_xscale("log")
        ax.set_xticks([0.32, 1, 10, 100, 1000])
        ax.set_xticklabels(["0", "1", "10", "100", "1000"])
        ax.axvline(0.62, color="0.75", lw=0.8, ls=":")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlabel("distance to nearest hotspot (kb)")

    # figure 1: pi vs distance
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for pop in POP_ORDER:
        t = binned[binned["pop"] == pop].sort_values("bin_idx")
        x = xs[t.bin_idx.to_numpy()]
        yerr = np.vstack([t.pi - t.ci_low, t.ci_high - t.pi])
        ax.errorbar(x, t.pi, yerr=yerr, marker="o", ms=4, lw=1.5, capsize=2,
                    color=POP_COLORS[pop], label=pop)
    style(ax)
    ax.set_ylabel(r"$\pi$")
    ax.set_title("Diversity vs distance to nearest recombination hotspot")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.png", dpi=200)
    plt.close(fig)

    # figure 2: the confound
    ref_bins = binned[binned["pop"] == POP_ORDER[0]].sort_values("bin_idx")
    x = xs[ref_bins.bin_idx.to_numpy()]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.9))
    axes[0].plot(x, ref_bins.median_local_rate, marker="o", ms=4, lw=1.5, color="#555555")
    axes[0].set_ylabel("local rate (cM/Mb), bin median")
    axes[0].set_title("Background recombination rate")
    axes[1].plot(x, ref_bins.median_callable, marker="o", ms=4, lw=1.5, color="#555555")
    axes[1].set_ylabel("callable sites per 2 kb window, median")
    axes[1].set_title("Callable-site density")
    for a in axes:
        style(a)
    fig.suptitle("What else changes with distance to a hotspot", y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_confounds.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # figure 3: stratified by local-rate quartile
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.9), sharey=True)
    qlab = ["Q1 (lowest rate)", "Q2", "Q3", "Q4 (highest rate)"]
    for q in range(4):
        ax = axes[q]
        for pop in POP_ORDER:
            t = strat[(strat.rate_q == q) & (strat["pop"] == pop)].sort_values("bin_idx")
            if t.empty:
                continue
            xq = xs[t.bin_idx.to_numpy()]
            yerr = np.vstack([t.pi - t.ci_low, t.ci_high - t.pi])
            ax.errorbar(xq, t.pi, yerr=yerr, marker="o", ms=3.5, lw=1.3, capsize=2,
                        color=POP_COLORS[pop], label=pop)
        style(ax)
        med_rate = strat[strat.rate_q == q].median_local_rate.median()
        ax.set_title(f"{qlab[q]}\nmedian {med_rate:.2f} cM/Mb", fontsize=10)
        if q == 0:
            ax.set_ylabel(r"$\pi$")
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Diversity vs hotspot distance within strata of local recombination rate",
                 y=1.03, fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_stratified.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"\nwrote {args.out_prefix}.png, _confounds.png, _stratified.png, "
          f"_bins.tsv, _stratified.tsv, _summary.tsv")


if __name__ == "__main__":
    main()
