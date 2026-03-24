#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_INPUT_PATH = Path("results/co_events_long.tsv")
DEFAULT_OUTPUT_PATH = Path("results/marey_map_co_events.png")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a Marey-style map from CO breakpoint calls.")
    parser.add_argument("input_path", nargs="?", default=str(DEFAULT_INPUT_PATH))
    parser.add_argument("output_path", nargs="?", default=str(DEFAULT_OUTPUT_PATH))
    return parser.parse_args()


def load_events(input_path: Path):
    per_chrom = defaultdict(list)
    chromosome_samples = defaultdict(set)

    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = int(row["chromosome"])
            bp = int(row["event_bp"])
            sample_id = row["sample_id"]
            per_chrom[chrom].append(bp)
            chromosome_samples[chrom].add(sample_id)

    return per_chrom, chromosome_samples


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_path)
    output_path = Path(args.output_path)
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/jri-arg-mpl")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    per_chrom, chromosome_samples = load_events(input_path)
    chromosomes = sorted(per_chrom)
    ncols = 2
    nrows = (len(chromosomes) + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 3.6 * nrows), constrained_layout=True)
    axes = axes.flatten()

    for ax, chrom in zip(axes, chromosomes):
        positions = sorted(per_chrom[chrom])
        n_samples = len(chromosome_samples[chrom])
        cumulative = []
        running = 0.0

        for _ in positions:
            running += 1.0 / n_samples
            cumulative.append(running)

        ax.step(positions, cumulative, where="post", color="#1f4e79", linewidth=1.4)
        ax.set_title(f"Chr {chrom}")
        ax.set_xlabel("Physical position (bp)")
        ax.set_ylabel("Cumulative COs per individual")
        ax.grid(alpha=0.25, linewidth=0.5)
        ax.ticklabel_format(style="sci", axis="x", scilimits=(6, 6))

    for ax in axes[len(chromosomes):]:
        ax.axis("off")

    fig.suptitle("Empirical Marey-Style Map from Called CO Breakpoints", fontsize=14)
    fig.savefig(output_path, dpi=200)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
