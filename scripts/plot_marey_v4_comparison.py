#!/usr/bin/env python3
"""Marey map comparing LR13/LR14 (maize) vs teosinte crossovers in AGPv4 coordinates."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent
LR_FILE = ROOT / "data/xo_Combined_LR13_LR14_parents2_AGPv4_filter0219.txt"
TEO_FILE = ROOT / "data/xo_ZeaGBSv27raw_RareAllelesC2TeoCurated_depth_AGPv4_filtered0210.txt"
OUT = ROOT / "results/marey_v4_maize_vs_teo.png"

# AGPv4 chromosome lengths (from chain file header / MaizeGDB)
CHR_LENGTHS = {
    1: 307_041_717,
    2: 244_442_276,
    3: 235_667_834,
    4: 246_994_605,
    5: 223_902_240,
    6: 174_033_170,
    7: 182_381_542,
    8: 181_122_637,
    9: 159_769_782,
    10: 150_982_314,
}
CHROMS = list(range(1, 11))


def load_xo(path, label):
    df = pd.read_csv(path, sep="\t")
    df["mid"] = (df["start"] + df["end"]) // 2
    df["label"] = label
    return df


def cumulative_co(df, chrom, n_indiv):
    """Return (positions, cumulative_COs_per_individual) for one chromosome."""
    sub = df[df["chr"] == chrom].copy()
    if sub.empty:
        return np.array([]), np.array([])
    mids = np.sort(sub["mid"].values)
    cumco = np.arange(1, len(mids) + 1) / n_indiv
    # Prepend 0 at position 0 and append final value at chr end
    pos = np.concatenate([[0], mids, [CHR_LENGTHS[chrom]]])
    vals = np.concatenate([[0], cumco, [cumco[-1]]])
    return pos, vals


def main():
    lr = load_xo(LR_FILE, "LR13/LR14 (maize)")
    teo = load_xo(TEO_FILE, "Teosinte")

    n_lr = lr["taxon"].nunique()
    n_teo = teo["taxon"].nunique()
    print(f"Maize taxons: {n_lr}, Teosinte taxons: {n_teo}")

    fig, axes = plt.subplots(2, 5, figsize=(18, 7), sharey=False)
    axes = axes.flatten()

    for ax, chrom in zip(axes, CHROMS):
        lr_pos, lr_vals = cumulative_co(lr, chrom, n_lr)
        teo_pos, teo_vals = cumulative_co(teo, chrom, n_teo)

        if len(lr_pos):
            ax.plot(lr_pos / 1e6, lr_vals, color="#2166ac", linewidth=1.0,
                    label=f"LR13/LR14 (maize, n={n_lr})")
        if len(teo_pos):
            ax.plot(teo_pos / 1e6, teo_vals, color="#d62728", linewidth=1.0,
                    label=f"Teosinte (n={n_teo})")

        ax.set_title(f"Chr{chrom}", fontsize=9)
        ax.set_xlabel("Position (Mb)", fontsize=7)
        ax.set_ylabel("Cumulative CO / individual", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.set_xlim(0, CHR_LENGTHS[chrom] / 1e6)

    import matplotlib.lines as mlines
    legend_handles = [
        mlines.Line2D([], [], color="#2166ac", linewidth=1.5,
                      label=f"LR13/LR14 maize (n={n_lr})"),
        mlines.Line2D([], [], color="#d62728", linewidth=1.5,
                      label=f"Teosinte (n={n_teo})"),
    ]
    fig.suptitle("Marey map: LR13/LR14 (maize) vs teosinte — AGPv4 coordinates", fontsize=11)
    fig.tight_layout()
    fig.legend(handles=legend_handles, loc="lower center", ncol=2,
               fontsize=9, frameon=True, framealpha=0.9,
               bbox_to_anchor=(0.5, -0.04))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
