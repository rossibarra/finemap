#!/usr/bin/env python3
"""Correlate 0-fold and 4-fold degenerate diversity with the FineMap recombination rate.

Same estimator, windows, bootstrap and partial correlation as
scripts/pi_vs_recombination.py -- this version reads two pi tables (one per
degeneracy class, from scripts/haploid_pi.py --restrict) and additionally
reports pi_0D / pi_4D, the ratio that indexes the efficacy of purifying
selection.  Under linked selection that ratio is expected to fall as
recombination rises, because selection is more efficient where linkage is
broken up faster.

Because 0D and 4D sites live only in coding sequence, a 100 kb window holds far
fewer callable sites than in the all-sites analysis, so the minimum-callable-site
filter has to be correspondingly lower (and is set separately per class).
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pi_vs_recombination import (  # noqa: E402
    BED_COLUMNS, POP_COLORS, spearman, partial_spearman, block_bootstrap_ci,
    gene_density, centromere_distance,
)

CLASSES = ["0D", "4D"]


def load_pi(path, bed, min_sites):
    pi = pd.read_csv(path, sep="\t")
    pi["start"] = pi.window_pos_1 - 1
    pi = pi.rename(columns={"chromosome": "chrom"})
    df = pi.merge(bed, on=["chrom", "start"], how="inner", validate="many_to_one")
    if df.empty:
        sys.exit(f"{path}: no windows matched the bed -- check chromosome naming")
    df = df[df.avg_pi.notna() & (df.no_sites >= min_sites) & df.rate.notna()].copy()
    return df.sort_values(["chrom", "start"]).reset_index(drop=True)


def ratio_se(b, n_boot=200, seed=1):
    """Bootstrap SD of the decile-aggregated pi_0D/pi_4D, resampling windows."""
    rng = np.random.default_rng(seed)
    d0 = b.count_diffs_0d.to_numpy(float); p0 = b.count_comparisons_0d.to_numpy(float)
    d4 = b.count_diffs_4d.to_numpy(float); p4 = b.count_comparisons_4d.to_numpy(float)
    vals = []
    for _ in range(n_boot):
        i = rng.integers(0, len(b), size=len(b))
        num, den = d4[i].sum(), p4[i].sum()
        if num <= 0 or den <= 0 or p0[i].sum() <= 0:
            continue
        vals.append((d0[i].sum() / p0[i].sum()) / (num / den))
    return float(np.std(vals)) if vals else np.nan


def decile_table(g, value_fn):
    q = pd.qcut(g.rate.rank(method="first"), 10, labels=False)
    rows = []
    for d, b in g.groupby(q):
        rows.append({"decile": int(d) + 1, "rate": b.rate.median(),
                     "n_windows": len(b), **value_fn(b)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pi-0d", required=True)
    ap.add_argument("--pi-4d", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--gff")
    ap.add_argument("--centromeres")
    ap.add_argument("--min-sites-0d", type=int, default=200)
    ap.add_argument("--min-sites-4d", type=int, default=50)
    ap.add_argument("--block-mb", type=float, default=5.0)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    bed = pd.read_csv(args.bed, sep="\t", header=None, names=BED_COLUMNS)
    bed["chrom"] = bed["chrom"].str.replace("Chr", "chr", regex=False)
    window_size = int((bed.end - bed.start).mode().iloc[0])
    block_windows = max(1, int(args.block_mb * 1e6 / window_size))

    tables = {
        "0D": load_pi(args.pi_0d, bed, args.min_sites_0d),
        "4D": load_pi(args.pi_4d, bed, args.min_sites_4d),
    }

    covariates = {}
    for cls, df in tables.items():
        names = []
        if args.gff:
            df["n_genes"] = gene_density(args.gff, df.chrom.to_numpy(),
                                         df.start.to_numpy(), window_size)
            names.append("n_genes")
        if args.centromeres:
            df["cen_dist"] = centromere_distance(args.centromeres, df.chrom.to_numpy(),
                                                 df.start.to_numpy(), df.end.to_numpy())
            if df.cen_dist.notna().all():
                names.append("cen_dist")
            else:
                df.drop(columns=["cen_dist"], inplace=True)
        covariates[cls] = names

    rows = []
    for cls, df in tables.items():
        for pop, g in df.groupby("pop"):
            rho = spearman(g.rate, g.avg_pi)
            lo, hi = block_bootstrap_ci(g, "rate", "avg_pi", block_windows,
                                        args.n_boot, args.seed)
            partial = (partial_spearman(g.rate, g.avg_pi,
                                        [g[c] for c in covariates[cls]])
                       if covariates[cls] else np.nan)
            rows.append({"class": cls, "pop": pop, "n_windows": len(g),
                         "spearman_rho": rho, "ci_low": lo, "ci_high": hi,
                         "partial_rho": partial,
                         "rho_rate_vs_callable_sites": spearman(g.rate, g.no_sites),
                         "mean_pi": g.count_diffs.sum() / g.count_comparisons.sum()})

    # --- paired windows for the 0D/4D ratio ---
    keys = ["pop", "chrom", "start", "end", "rate"]
    merged = tables["0D"][keys + ["avg_pi", "count_diffs", "count_comparisons", "no_sites"]].merge(
        tables["4D"][keys + ["avg_pi", "count_diffs", "count_comparisons", "no_sites"]],
        on=keys, suffixes=("_0d", "_4d"), validate="one_to_one")
    merged = merged[merged.avg_pi_4d > 0].copy()
    merged["ratio"] = merged.avg_pi_0d / merged.avg_pi_4d
    if args.gff:
        merged["n_genes"] = gene_density(args.gff, merged.chrom.to_numpy(),
                                         merged.start.to_numpy(), window_size)
    if args.centromeres:
        merged["cen_dist"] = centromere_distance(args.centromeres, merged.chrom.to_numpy(),
                                                 merged.start.to_numpy(), merged.end.to_numpy())
    ratio_cov = [c for c in ("n_genes", "cen_dist") if c in merged.columns]
    merged = merged.sort_values(["chrom", "start"]).reset_index(drop=True)

    for pop, g in merged.groupby("pop"):
        rho = spearman(g.rate, g.ratio)
        lo, hi = block_bootstrap_ci(g, "rate", "ratio", block_windows,
                                    args.n_boot, args.seed)
        partial = (partial_spearman(g.rate, g.ratio, [g[c] for c in ratio_cov])
                   if ratio_cov else np.nan)
        rows.append({"class": "0D/4D", "pop": pop, "n_windows": len(g),
                     "spearman_rho": rho, "ci_low": lo, "ci_high": hi,
                     "partial_rho": partial, "rho_rate_vs_callable_sites": np.nan,
                     "mean_pi": ((g.count_diffs_0d.sum() / g.count_comparisons_0d.sum()) /
                                 (g.count_diffs_4d.sum() / g.count_comparisons_4d.sum()))})

    summary = pd.DataFrame(rows)
    summary.to_csv(f"{args.out_prefix}_summary.tsv", sep="\t", index=False,
                   float_format="%.4f")

    # --- decile tables ---
    dec_rows = []
    for cls, df in tables.items():
        for pop, g in df.groupby("pop"):
            t = decile_table(g, lambda b: {
                "pi": b.count_diffs.sum() / b.count_comparisons.sum(),
                "se": b.avg_pi.std() / np.sqrt(len(b)),
                "median_sites": b.no_sites.median()})
            t.insert(0, "pop", pop)
            t.insert(0, "class", cls)
            dec_rows.append(t)
    for pop, g in merged.groupby("pop"):
        t = decile_table(g, lambda b: {
            "pi": ((b.count_diffs_0d.sum() / b.count_comparisons_0d.sum()) /
                   (b.count_diffs_4d.sum() / b.count_comparisons_4d.sum())),
            "se": ratio_se(b, seed=args.seed),
            "median_sites": b.no_sites_4d.median()})
        t.insert(0, "pop", pop)
        t.insert(0, "class", "0D/4D")
        dec_rows.append(t)
    deciles = pd.concat(dec_rows, ignore_index=True)
    deciles.to_csv(f"{args.out_prefix}_deciles.tsv", sep="\t", index=False,
                   float_format="%.6g")
    merged.to_csv(f"{args.out_prefix}_windows.tsv.gz", sep="\t", index=False)

    for cls in CLASSES:
        n = len(tables[cls]) // max(tables[cls]["pop"].nunique(), 1)
        thresh = args.min_sites_0d if cls == "0D" else args.min_sites_4d
        print(f"{cls}: {n} of {len(bed)} windows pass >= {thresh} callable sites")
    print(f"paired (both classes pass, pi_4D > 0): "
          f"{len(merged) // max(merged['pop'].nunique(), 1)} windows")
    print()
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    pops = sorted(tables["0D"]["pop"].unique())

    # --- figure 1: pi by recombination decile, 0D and 4D ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, cls in zip(axes, CLASSES):
        for pop in pops:
            t = deciles[(deciles["class"] == cls) & (deciles["pop"] == pop)]
            ax.errorbar(t.rate, t.pi, yerr=t.se, marker="o", ms=4, lw=1.5,
                        capsize=2, label=pop, color=POP_COLORS.get(pop))
        ax.set_xlabel("recombination rate (cM/Mb), decile median")
        ax.set_ylabel(rf"$\pi_{{{cls}}}$")
        label = "0-fold degenerate (nonsynonymous proxy)" if cls == "0D" \
            else "4-fold degenerate (synonymous proxy)"
        ax.set_title(label, fontsize=11)
    axes[0].legend(frameon=False, fontsize=9)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.png", dpi=200)

    # --- figure 2: the 0D/4D ratio ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    for pop in pops:
        t = deciles[(deciles["class"] == "0D/4D") & (deciles["pop"] == pop)]
        ax.errorbar(t.rate, t.pi, yerr=t.se, marker="o", ms=4, lw=1.5, capsize=2,
                    label=pop, color=POP_COLORS.get(pop))
    ax.set_xlabel("recombination rate (cM/Mb), decile median")
    ax.set_ylabel(r"$\pi_{0D} / \pi_{4D}$")
    ax.set_title("Efficacy of purifying selection by decile", fontsize=11)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    for pop in pops:
        g = merged[merged["pop"] == pop]
        ax.scatter(g.rate, g.ratio, s=2, alpha=0.12, color=POP_COLORS.get(pop),
                   label=pop, rasterized=True)
    ax.set_xlabel("recombination rate (cM/Mb)")
    ax.set_ylabel(r"$\pi_{0D} / \pi_{4D}$")
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set_ylim(0, np.nanquantile(merged.ratio, 0.995))
    ax.set_title(f"All paired windows ({window_size // 1000} kb)", fontsize=11)
    leg = ax.legend(frameon=False, fontsize=9, markerscale=6)
    for h in leg.legend_handles:
        h.set_alpha(1)

    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_ratio.png", dpi=200)

    print(f"\nwrote {args.out_prefix}.png, {args.out_prefix}_ratio.png, "
          f"{args.out_prefix}_summary.tsv, {args.out_prefix}_deciles.tsv, "
          f"{args.out_prefix}_windows.tsv.gz")


if __name__ == "__main__":
    main()
