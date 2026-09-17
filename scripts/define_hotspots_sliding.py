#!/usr/bin/env python3
"""Define recombination hotspots from a sliding-window average of the FineMap map.

The interval-density map (data/finemap_v5.bed) is piecewise constant over 214k
variable-width intervals, many only a few bp wide.  Thresholding those intervals
directly rewards narrowness: rate is cM per bp, so a crossover localized to 40 bp
scores enormously whether or not the region is genuinely hot.

Averaging the rate over a fixed 1 kb window removes that sensitivity -- a single
narrow spike is diluted by its neighbours, and only sustained elevation survives.
The window average is computed exactly rather than by binning: the map defines a
piecewise-linear cumulative genetic position, so the mean rate over [a, b] is

    (cM(b) - cM(a)) / ((b - a) / 1e6)

with cM() linearly interpolated at the interval breakpoints.  Windows are stepped
along each chromosome, those above the threshold are kept, and overlapping or
contiguous passing windows are merged into hotspot regions.
"""

import argparse

import numpy as np
import pandas as pd

COLS = ["chrom", "start", "end", "cm_start", "cm_end", "rate"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, help="FineMap interval-density BED")
    ap.add_argument("--window", type=int, default=1000, help="sliding window width in bp")
    ap.add_argument("--step", type=int, default=100, help="slide step in bp")
    ap.add_argument("--fold", type=float, default=20.0,
                    help="threshold as a multiple of the length-weighted genome mean")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    bed = pd.read_csv(args.bed, sep="\t", header=None, names=COLS)
    bed["chrom"] = bed.chrom.str.replace("Chr", "chr", regex=False)
    bed = bed.sort_values(["chrom", "start"])

    total_bp = int((bed.end - bed.start).sum())
    total_cm = float((bed.cm_end - bed.cm_start).sum())
    mean_rate = total_cm / (total_bp / 1e6)
    thr = args.fold * mean_rate
    print(f"length-weighted genome mean  {mean_rate:.4f} cM/Mb")
    print(f"threshold ({args.fold:g}x)            {thr:.4f} cM/Mb")
    print(f"window / step                {args.window} / {args.step} bp\n")

    regions = []
    n_windows_total = 0
    n_pass_total = 0

    for chrom, g in bed.groupby("chrom", sort=False):
        # Cumulative genetic position at each interval breakpoint. cm_start/cm_end are
        # already cumulative within a chromosome, so the breakpoints and their cM values
        # define the piecewise-linear map directly.
        pos = np.concatenate([g.start.to_numpy(), g.end.to_numpy()[-1:]]).astype(np.float64)
        cm = np.concatenate([g.cm_start.to_numpy(), g.cm_end.to_numpy()[-1:]]).astype(np.float64)

        first, last = pos[0], pos[-1]
        if last - first < args.window:
            continue
        starts = np.arange(first, last - args.window + 1, args.step, dtype=np.float64)
        ends = starts + args.window
        rate = (np.interp(ends, pos, cm) - np.interp(starts, pos, cm)) / (args.window / 1e6)

        n_windows_total += len(starts)
        hot = rate > thr
        n_pass_total += int(hot.sum())
        if not hot.any():
            continue

        # Merge passing windows that overlap or abut.
        hs = starts[hot].astype(np.int64)
        he = ends[hot].astype(np.int64)
        breaks = np.flatnonzero(hs[1:] > he[:-1])
        grp_start = np.concatenate([[0], breaks + 1])
        grp_end = np.concatenate([breaks, [len(hs) - 1]])
        for a, b in zip(grp_start, grp_end):
            s, e = int(hs[a]), int(he[b])
            peak = float(rate[hot][a:b + 1].max())
            mean_cm = float(np.interp(e, pos, cm) - np.interp(s, pos, cm))
            regions.append((chrom, s, e, e - s, mean_cm, mean_cm / ((e - s) / 1e6), peak))

    out = pd.DataFrame(regions, columns=["chrom", "start", "end", "width",
                                         "cm", "mean_rate", "peak_window_rate"])
    out["chrom_n"] = out.chrom.str.replace("chr", "", regex=False).astype(int)
    out = out.sort_values(["chrom_n", "start"]).drop(columns="chrom_n")

    print(f"windows evaluated            {n_windows_total:,}")
    print(f"windows above threshold      {n_pass_total:,} "
          f"({100 * n_pass_total / n_windows_total:.3f}%)")
    print(f"merged hotspot regions       {len(out):,}\n")
    print("width (bp):")
    print(out.width.describe(percentiles=[.05, .25, .5, .75, .95]).round(0).to_string())
    print(f"\ntotal hotspot bp             {int(out.width.sum()):,} "
          f"({100 * out.width.sum() / total_bp:.4f}% of mapped sequence)")
    print(f"total hotspot cM             {out.cm.sum():.2f} "
          f"({100 * out.cm.sum() / total_cm:.2f}% of the map)")
    print(f"cM enrichment                {(out.cm.sum() / total_cm) / (out.width.sum() / total_bp):.1f}x")
    print("\nregions per chromosome:")
    print(out.groupby("chrom").size().to_string())

    out[["chrom", "start", "end", "mean_rate"]].to_csv(
        args.out, sep="\t", header=False, index=False, float_format="%.4f")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
