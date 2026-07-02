#!/usr/bin/env python3
"""Plot recombination rate (cM/Mb) along each chromosome for finemap_v5 and Ogut."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent
FINEMAP = ROOT / "data/finemap_v5.bed"
OGUT_V5 = ROOT / "data/ogut_v5.csv"
CENTROMERES = ROOT / "data/NAM_centromere_coords-cenH3.csv"
OUT = ROOT / "results/rate_along_chromosomes.png"

CHROMS = [f"Chr{i}" for i in range(1, 11)]


def load_centromeres():
    cen = pd.read_csv(CENTROMERES)
    cen["Chr"] = cen["Chr"].str.replace(r"^chr", "Chr", regex=True)
    return cen.set_index("Chr")[["Start", "End"]]


def load_finemap():
    fm = pd.read_csv(
        FINEMAP, sep="\t", header=None,
        names=["chr", "start", "end", "cM_start", "cM_end", "cM_per_Mb"]
    )
    return fm


def load_ogut():
    og = pd.read_csv(OGUT_V5)
    og = og.rename(columns={"chr": "chrom", "pos_v5": "pos"})
    # Compute per-interval rate and midpoint from consecutive markers per chrom
    rows = []
    for chrom, grp in og.groupby("chrom"):
        grp = grp.sort_values("pos").reset_index(drop=True)
        mid = (grp["pos"].values[:-1] + grp["pos"].values[1:]) / 2
        dcm = np.diff(grp["cM_norm"].values)
        dbp = np.diff(grp["pos"].values)
        with np.errstate(divide="ignore", invalid="ignore"):
            rate = np.where(dbp > 0, dcm / (dbp / 1e6), np.nan)
        rows.append(pd.DataFrame({"chrom": chrom, "mid": mid, "rate": rate}))
    return pd.concat(rows, ignore_index=True)


def plot(finemap, ogut, centromeres):
    fig, axes = plt.subplots(2, 5, figsize=(18, 7), sharey=False)
    axes = axes.flatten()

    for ax, chrom in zip(axes, CHROMS):
        num = int(chrom.replace("Chr", ""))

        fm = finemap[finemap["chr"] == chrom].sort_values("start")
        og = ogut[ogut["chrom"] == chrom].sort_values("mid")
        og_plot = og[og["rate"] >= 0]  # drop negative rates from liftover artefacts

        if fm.empty:
            ax.set_visible(False)
            continue

        # Ogut: scatter of interval rates
        ax.scatter(
            og_plot["mid"] / 1e6, og_plot["rate"],
            s=1, color="#999999", alpha=0.6, zorder=2, label="Ogut (lifted)"
        )

        # finemap_v5: step line using segment midpoints and cM_per_Mb
        mid_fm = (fm["start"].values + fm["end"].values) / 2
        ax.step(
            mid_fm / 1e6, fm["cM_per_Mb"].values,
            where="mid", color="#d62728", linewidth=0.8, zorder=3, label="finemap_v5"
        )

        if chrom in centromeres.index:
            ax.axvspan(
                centromeres.loc[chrom, "Start"] / 1e6,
                centromeres.loc[chrom, "End"] / 1e6,
                color="#aec6e8", alpha=0.5, zorder=1
            )

        ax.set_title(f"Chr{num}", fontsize=9)
        ax.set_xlabel("Position (Mb)", fontsize=7)
        ax.set_ylabel("cM/Mb", fontsize=7)
        ax.tick_params(labelsize=6)

    legend_handles = [
        mlines.Line2D([], [], color="#999999", marker="o", markersize=4,
                      linestyle="None", label="Ogut (v5 lifted)"),
        mlines.Line2D([], [], color="#d62728", linewidth=1.5, label="finemap_v5"),
        mpatches.Patch(facecolor="#aec6e8", alpha=0.5, label="Centromere (B73)"),
    ]
    fig.suptitle("Recombination rate (cM/Mb) along chromosomes", fontsize=11)
    fig.tight_layout()
    fig.legend(handles=legend_handles, loc="lower center", ncol=3,
               fontsize=9, frameon=True, framealpha=0.9,
               bbox_to_anchor=(0.5, -0.03))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"Wrote {OUT}")


def main():
    finemap = load_finemap()
    ogut = load_ogut()
    centromeres = load_centromeres()
    plot(finemap, ogut, centromeres)


if __name__ == "__main__":
    main()
