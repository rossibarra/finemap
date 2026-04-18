#!/usr/bin/env python3
"""Marey map comparing Ogut v5-lifted positions vs finemap_v5.bed."""

import subprocess
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent
OGUT_V2 = ROOT / "data/ogut_fifthcM_map_agpv2.csv"
CHAIN = ROOT / "data/v2v5.chain"
FINEMAP = ROOT / "data/finemap_v5.bed"
OUT = ROOT / "results/marey_ogut_vs_finemap.png"

CHROMS = [f"Chr{i}" for i in range(1, 11)]
CENTROMERES = ROOT / "data/NAM_centromere_coords-cenH3.csv"


def load_centromeres():
    cen = pd.read_csv(CENTROMERES)
    cen["Chr"] = cen["Chr"].str.replace(r"^chr", "Chr", regex=True)
    return cen.set_index("Chr")[["Start", "End"]]


def lift_ogut():
    ogut = pd.read_csv(OGUT_V2)
    # CrossMap expects BED: chr start end name
    # Ogut positions are in AGPv2 (bare numbers); v2v5.chain also uses bare numbers
    ogut["chr_str"] = ogut["chromosome"].astype(str)
    ogut["end"] = ogut["position"] + 1

    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False) as fh:
        tmp_in = fh.name
        for _, row in ogut.iterrows():
            fh.write(f"{row.chr_str}\t{row.position}\t{row.end}\t{row.SNP_newID}\n")

    tmp_out = tmp_in.replace(".bed", "_v5.bed")
    crossmap = subprocess.run(["which", "CrossMap"], capture_output=True, text=True).stdout.strip()
    if not crossmap:
        crossmap = "/opt/anaconda3/bin/CrossMap"
    subprocess.run(
        [crossmap, "bed", str(CHAIN), tmp_in, tmp_out],
        check=True,
        capture_output=True,
    )

    lifted = pd.read_csv(
        tmp_out, sep="\t", header=None, names=["chr_v5", "start_v5", "end_v5", "SNP_newID"]
    )
    # chain already outputs Chr1..Chr10; keep as-is
    merged = ogut[["SNP_newID", "chromosome", "cM"]].merge(lifted, on="SNP_newID", how="inner")
    # Shift cM per chromosome so min = 0
    merged["cM_norm"] = merged.groupby("chromosome")["cM"].transform(lambda x: x - x.min())
    return merged


def load_finemap():
    fm = pd.read_csv(
        FINEMAP, sep="\t", header=None,
        names=["chr", "start", "end", "cM_start", "cM_end", "cM_per_Mb"]
    )
    return fm


def plot(ogut_v5, finemap, centromeres):
    fig, axes = plt.subplots(2, 5, figsize=(18, 7), sharey=False)
    axes = axes.flatten()

    for ax, chrom in zip(axes, CHROMS):
        num = int(chrom.replace("Chr", ""))

        og = ogut_v5[ogut_v5["chr_v5"] == chrom].sort_values("start_v5")
        fm = finemap[finemap["chr"] == chrom].sort_values("start")

        if og.empty or fm.empty:
            ax.set_visible(False)
            continue

        # Ogut: scatter of lifted marker positions
        ax.scatter(
            og["start_v5"] / 1e6, og["cM_norm"],
            s=2, color="#999999", zorder=2, label="Ogut (lifted)"
        )

        # finemap_v5: step line from segment breakpoints
        xs = []
        ys = []
        for _, row in fm.iterrows():
            xs += [row["start"] / 1e6, row["end"] / 1e6]
            ys += [row["cM_start"], row["cM_end"]]
        ax.plot(xs, ys, color="#d62728", linewidth=1.0, zorder=3, label="finemap_v5")

        if chrom in centromeres.index:
            cen_start = centromeres.loc[chrom, "Start"] / 1e6
            cen_end = centromeres.loc[chrom, "End"] / 1e6
            ax.axvspan(cen_start, cen_end, color="#aec6e8", alpha=0.5, zorder=1,
                       label="Centromere")

        ax.set_title(f"Chr{num}", fontsize=9)
        ax.set_xlabel("Position (Mb)", fontsize=7)
        ax.set_ylabel("cM", fontsize=7)
        ax.tick_params(labelsize=6)

    import matplotlib.lines as mlines
    import matplotlib.patches as mpatches
    legend_handles = [
        mlines.Line2D([], [], color="#999999", marker="o", markersize=4,
                      linestyle="None", label="Ogut (v5 lifted)"),
        mlines.Line2D([], [], color="#d62728", linewidth=1.5, label="finemap_v5"),
        mpatches.Patch(facecolor="#aec6e8", alpha=0.5, label="Centromere (B73)"),
    ]
    fig.suptitle("Marey map: Ogut (v5 lifted) vs finemap_v5", fontsize=11)
    fig.tight_layout()
    fig.legend(handles=legend_handles, loc="lower center", ncol=2,
               fontsize=9, frameon=True, framealpha=0.9,
               bbox_to_anchor=(0.5, -0.04))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"Wrote {OUT}")


def main():
    print("Lifting Ogut markers to v5...")
    ogut_v5 = lift_ogut()
    print(f"  {len(ogut_v5)} markers lifted")
    finemap = load_finemap()
    centromeres = load_centromeres()
    print("Plotting...")
    plot(ogut_v5, finemap, centromeres)


if __name__ == "__main__":
    main()
