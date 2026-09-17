#!/usr/bin/env python3
"""Define recombination hotspots from the FineMap interval-density map.

The interval map (data/finemap_v5.bed) is a piecewise-constant cM/Mb function over
214k variable-width intervals.  A hotspot is a run of contiguous intervals whose rate
exceeds a multiple of the *length-weighted* genome-wide mean (total cM / total Mb).
The unweighted mean of interval rates is badly inflated because the map is dominated
by very short intervals, which mechanically carry huge rates.

Merged regions narrower than --min-width are dropped: a crossover localized to a few
bp yields an enormous apparent rate that is a resolution artifact, not evidence of a
hotspot.
"""

import argparse
import numpy as np
import pandas as pd

COLS = ["chrom", "start", "end", "cm_start", "cm_end", "rate"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--fold", type=float, default=30.0)
    ap.add_argument("--min-width", type=int, default=50)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    bed = pd.read_csv(args.bed, sep="\t", header=None, names=COLS)
    bed["chrom"] = bed.chrom.str.replace("Chr", "chr", regex=False)
    bed["width"] = bed.end - bed.start
    bed["cm"] = bed.cm_end - bed.cm_start

    total_bp = int(bed.width.sum())
    total_cm = float(bed.cm.sum())
    weighted_mean = total_cm / (total_bp / 1e6)
    unweighted = float(bed.rate.mean())
    print(f"intervals                  {len(bed):,}")
    print(f"total bp                   {total_bp:,}")
    print(f"total cM                   {total_cm:.1f}")
    print(f"length-weighted mean rate  {weighted_mean:.4f} cM/Mb")
    print(f"unweighted mean of rates   {unweighted:.4f} cM/Mb")
    print(f"median interval width      {bed.width.median():.0f} bp")

    thresh = args.fold * weighted_mean
    print(f"\nthreshold ({args.fold:g}x)          {thresh:.3f} cM/Mb")
    hot = bed[bed.rate >= thresh].copy()
    print(f"intervals above threshold  {len(hot):,}")
    print(f"  median width             {hot.width.median():.0f} bp")
    print(f"  fraction <= 100 bp       {(hot.width <= 100).mean():.3f}")

    # merge contiguous runs within a chromosome
    hot = hot.sort_values(["chrom", "start"]).reset_index(drop=True)
    new_group = (hot.chrom != hot.chrom.shift()) | (hot.start != hot.end.shift())
    hot["grp"] = new_group.cumsum()
    merged = hot.groupby("grp").agg(chrom=("chrom", "first"), start=("start", "min"),
                                    end=("end", "max"), cm=("cm", "sum"),
                                    n_intervals=("rate", "size")).reset_index(drop=True)
    merged["width"] = merged.end - merged.start
    merged["rate"] = merged.cm / (merged.width / 1e6)
    print(f"merged regions             {len(merged):,}")

    keep = merged[merged.width >= args.min_width].copy()
    print(f"after >= {args.min_width} bp filter      {len(keep):,}")
    print(f"\nsurviving hotspot width distribution (bp):")
    for q in [0, 5, 25, 50, 75, 95, 100]:
        print(f"  p{q:<3d} {np.percentile(keep.width, q):>12,.0f}")
    print(f"  mean  {keep.width.mean():>12,.0f}")
    print(f"  fraction at exactly the {args.min_width} bp floor: "
          f"{(keep.width == args.min_width).mean():.3f}")
    print(f"  fraction <= 100 bp: {(keep.width <= 100).mean():.3f}")
    print(f"  fraction <= 1 kb:   {(keep.width <= 1000).mean():.3f}")
    print(f"\ntotal hotspot bp           {int(keep.width.sum()):,} "
          f"({100 * keep.width.sum() / total_bp:.3f}% of the map)")
    print(f"total hotspot cM           {keep.cm.sum():.1f} "
          f"({100 * keep.cm.sum() / total_cm:.2f}% of the map)")
    print(f"mean hotspot rate          {keep.cm.sum() / (keep.width.sum() / 1e6):.1f} cM/Mb")

    print("\nper chromosome:")
    per = keep.groupby("chrom").agg(n=("width", "size"), bp=("width", "sum"),
                                    cm=("cm", "sum"), median_width=("width", "median"))
    per = per.reindex([f"chr{i}" for i in range(1, 11)])
    print(per.to_string(float_format=lambda v: f"{v:,.1f}"))

    # artifact diagnostic: does rate simply track 1/width?
    def spearman(x, y):
        rx = np.array(pd.Series(x).rank(), dtype=float)
        ry = np.array(pd.Series(y).rank(), dtype=float)
        rx -= rx.mean()
        ry -= ry.mean()
        return float((rx * ry).sum() / np.sqrt((rx * rx).sum() * (ry * ry).sum()))

    print(f"\nrank corr(interval rate, 1/width), all {len(bed):,} intervals: "
          f"{spearman(bed.rate, 1.0 / bed.width):.3f}")
    sel = bed[bed.rate >= thresh]
    print(f"rank corr(interval rate, 1/width), above threshold:          "
          f"{spearman(sel.rate, 1.0 / sel.width):.3f}")
    print(f"rank corr(hotspot rate, 1/width), surviving hotspots:        "
          f"{spearman(keep.rate, 1.0 / keep.width):.3f}")

    keep[["chrom", "start", "end", "cm", "rate", "n_intervals"]].to_csv(
        args.out, sep="\t", header=False, index=False, float_format="%.6g")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
