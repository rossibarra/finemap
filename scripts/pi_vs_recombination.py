#!/usr/bin/env python3
"""Correlate windowed nucleotide diversity with the FineMap recombination rate.

Joins a pixy-format pi table to the piecewise-constant map and reports, per
population, the Spearman correlation between pi and cM/Mb.

Both variables are strongly autocorrelated along a chromosome, so the usual
asymptotic p-value over ~21k windows is meaningless -- adjacent windows are not
independent observations.  Uncertainty here comes from a moving-block bootstrap
that resamples contiguous blocks of windows, which preserves that autocorrelation.

Recombination, gene density and pi all covary with distance to the centromere,
so a partial correlation controlling for gene density and centromere distance is
reported alongside the raw one.
"""

import argparse
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BED_COLUMNS = ["chrom", "start", "end", "cm_start", "cm_end", "rate"]
POP_COLORS = {"maize": "#B8860B", "parviglumis": "#2E7D5B", "mexicana": "#3B6FA0"}


def spearman(x, y):
    """Spearman rho without a scipy dependency: Pearson on ranks."""
    rx = np.array(pd.Series(x).rank(), dtype=float)
    ry = np.array(pd.Series(y).rank(), dtype=float)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / denom) if denom > 0 else np.nan


def partial_spearman(x, y, covariates):
    """Spearman of x and y after linearly regressing rank-transformed covariates out."""
    rx = np.array(pd.Series(x).rank(), dtype=float)
    ry = np.array(pd.Series(y).rank(), dtype=float)
    cov = np.column_stack([np.array(pd.Series(c).rank(), dtype=float) for c in covariates])
    design = np.column_stack([np.ones(len(rx)), cov])
    rx_resid = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    ry_resid = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    return spearman(rx_resid, ry_resid)


def block_bootstrap_ci(df, xcol, ycol, block_windows, n_boot, seed, alpha=0.05):
    """Percentile CI for Spearman rho from a moving-block bootstrap within chromosomes."""
    rng = np.random.default_rng(seed)
    blocks = []
    for _, g in df.groupby("chrom", sort=False):
        x = np.array(g[xcol], dtype=float)
        y = np.array(g[ycol], dtype=float)
        for i in range(0, len(g), block_windows):
            if i + 2 <= len(g):
                blocks.append((x[i:i + block_windows], y[i:i + block_windows]))
    if len(blocks) < 2:
        return (np.nan, np.nan)

    n_target = len(df)
    rhos = []
    for _ in range(n_boot):
        picks = rng.integers(0, len(blocks), size=int(np.ceil(n_target / block_windows)))
        xs = np.concatenate([blocks[p][0] for p in picks])
        ys = np.concatenate([blocks[p][1] for p in picks])
        rho = spearman(xs, ys)
        if not np.isnan(rho):
            rhos.append(rho)
    if not rhos:
        return (np.nan, np.nan)
    return (float(np.quantile(rhos, alpha / 2)), float(np.quantile(rhos, 1 - alpha / 2)))


def gene_density(gff_path, chrom, start, window_size):
    """Genes per window, counted by gene midpoint -- the same convention as
    scripts/plot_rate_vs_gene_density.py, so the two analyses stay comparable."""
    counts = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            key = f[0].lower()
            mid = (int(f[3]) + int(f[4])) // 2
            counts[(key, (mid - 1) // window_size)] = counts.get((key, (mid - 1) // window_size), 0) + 1
    return np.array([counts.get((c.lower(), s // window_size), 0)
                     for c, s in zip(chrom, start)], dtype=float)


def centromere_distance(csv_path, chrom, start, end):
    """Absolute bp from window midpoint to the centromere midpoint."""
    cen = pd.read_csv(csv_path)
    cols = {c.lower(): c for c in cen.columns}
    chrom_col = next((cols[c] for c in cols if c in ("chr", "chrom", "chromosome")), None)
    start_col = next((cols[c] for c in cols if "start" in c), None)
    end_col = next((cols[c] for c in cols if "end" in c or "stop" in c), None)
    if not all([chrom_col, start_col, end_col]):
        sys.exit(f"could not identify chrom/start/end columns in {csv_path}: {list(cen.columns)}")
    mid = {}
    for _, r in cen.iterrows():
        key = str(r[chrom_col])
        key = key if key.lower().startswith("chr") else f"chr{key}"
        key = key.replace("Chr", "chr")
        mid.setdefault(key, []).append((float(r[start_col]) + float(r[end_col])) / 2)
    centers = {k: float(np.mean(v)) for k, v in mid.items()}
    win_mid = (start + end) / 2
    return np.array([abs(m - centers.get(c, np.nan)) for c, m in zip(chrom, win_mid)])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pi", required=True, help="pixy-format pi table")
    ap.add_argument("--bed", required=True, help="FineMap bed (chrom start end cM_start cM_end cM/Mb)")
    ap.add_argument("--gff", help="gene annotation for the partial correlation")
    ap.add_argument("--centromeres", help="centromere CSV for the partial correlation")
    ap.add_argument("--min-sites", type=int, default=5000,
                    help="drop windows with fewer callable sites (default 5000)")
    ap.add_argument("--block-mb", type=float, default=5.0, help="bootstrap block size in Mb")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    bed = pd.read_csv(args.bed, sep="\t", header=None, names=BED_COLUMNS)
    bed["chrom"] = bed["chrom"].str.replace("Chr", "chr", regex=False)
    window_size = int((bed.end - bed.start).mode().iloc[0])

    pi = pd.read_csv(args.pi, sep="\t")
    pi["start"] = pi.window_pos_1 - 1
    pi = pi.rename(columns={"chromosome": "chrom"})

    df = pi.merge(bed, on=["chrom", "start"], how="inner", validate="many_to_one")
    if df.empty:
        sys.exit("no windows matched between the pi table and the bed -- check chromosome naming")

    df = df[df.avg_pi.notna() & (df.no_sites >= args.min_sites) & df.rate.notna()].copy()

    covariate_names = []
    if args.gff:
        df["n_genes"] = gene_density(args.gff, df.chrom.to_numpy(),
                                     df.start.to_numpy(), window_size)
        covariate_names.append("n_genes")
    if args.centromeres:
        df["cen_dist"] = centromere_distance(args.centromeres, df.chrom.to_numpy(),
                                             df.start.to_numpy(), df.end.to_numpy())
        if df.cen_dist.notna().all():
            covariate_names.append("cen_dist")
        else:
            print("warning: centromere coords missing for some chromosomes; skipping that covariate",
                  file=sys.stderr)
            df = df.drop(columns=["cen_dist"])

    df = df.sort_values(["chrom", "start"]).reset_index(drop=True)
    block_windows = max(1, int(args.block_mb * 1e6 / window_size))

    rows = []
    for pop, g in df.groupby("pop"):
        rho = spearman(g.rate, g.avg_pi)
        lo, hi = block_bootstrap_ci(g, "rate", "avg_pi", block_windows, args.n_boot, args.seed)
        partial = (partial_spearman(g.rate, g.avg_pi, [g[c] for c in covariate_names])
                   if covariate_names else np.nan)
        rho_sites = spearman(g.rate, g.no_sites)
        rows.append({"pop": pop, "n_windows": len(g), "spearman_rho": rho,
                     "ci_low": lo, "ci_high": hi, "partial_rho": partial,
                     "rho_rate_vs_callable_sites": rho_sites,
                     "mean_pi": g.count_diffs.sum() / g.count_comparisons.sum()})

    summary = pd.DataFrame(rows).sort_values("pop")
    summary.to_csv(f"{args.out_prefix}_summary.tsv", sep="\t", index=False, float_format="%.4f")
    df.to_csv(f"{args.out_prefix}_windows.tsv.gz", sep="\t", index=False)

    print(f"\nwindows retained (>= {args.min_sites} callable sites): "
          f"{len(df)//df['pop'].nunique()} of {len(bed)}")
    if covariate_names:
        print(f"partial correlation controls for: {', '.join(covariate_names)}")
    print()
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    # --- figure: pi by recombination decile, plus the raw cloud ---
    pops = list(summary["pop"])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    for pop in pops:
        g = df[df["pop"] == pop]
        q = pd.qcut(g.rate.rank(method="first"), 10, labels=False)
        binned = g.groupby(q).apply(
            lambda b: pd.Series({
                "rate": b.rate.median(),
                "pi": b.count_diffs.sum() / b.count_comparisons.sum(),
                "se": b.avg_pi.std() / np.sqrt(len(b)),
            }), include_groups=False)
        ax.errorbar(binned.rate, binned.pi, yerr=binned.se, marker="o", ms=4,
                    lw=1.5, capsize=2, label=pop, color=POP_COLORS.get(pop))
    ax.set_xlabel("recombination rate (cM/Mb), decile median")
    ax.set_ylabel(r"$\pi$")
    ax.set_title("Diversity by recombination decile")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    for pop in pops:
        g = df[df["pop"] == pop]
        ax.scatter(g.rate, g.avg_pi, s=2, alpha=0.12, color=POP_COLORS.get(pop),
                   label=pop, rasterized=True)
    ax.set_xlabel("recombination rate (cM/Mb)")
    ax.set_ylabel(r"$\pi$")
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set_title(f"All windows ({window_size // 1000} kb)")
    leg = ax.legend(frameon=False, fontsize=9, markerscale=6)
    for h in leg.legend_handles:
        h.set_alpha(1)

    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.png", dpi=200)
    print(f"\nwrote {args.out_prefix}.png, {args.out_prefix}_summary.tsv, "
          f"{args.out_prefix}_windows.tsv.gz")


if __name__ == "__main__":
    main()
