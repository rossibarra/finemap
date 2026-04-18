#!/usr/bin/env python3
"""Run metaplot on the four example BED files and combine into a 2×2 figure."""

import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent.parent

EXAMPLES = [
    ("data/example1.bed", "Uniform genome-wide"),
    ("data/example2.bed", "5′ gene flanks"),
    ("data/example3.bed", "Gene bodies"),
    ("data/example4.bed", "3′ gene flanks"),
]

OUT = ROOT / "results/example_metaplots.png"


def run_metaplot(bed_path, title, out_path):
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/metaplot.py"),
            "--gff", str(ROOT / "data/v5.genes.gff3"),
            "--input", str(bed_path),
            "--bin-size", "100",
            "--flanking-bp", "5000",
            "--body-bins", "5",
            "--title", title,
            "--output", str(out_path),
            "--no-show",
        ],
        check=True,
    )


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        panels = []
        for i, (bed, title) in enumerate(EXAMPLES):
            out_png = tmpdir / f"panel{i+1}.png"
            print(f"  [{i+1}/4] {title} ({bed})")
            run_metaplot(ROOT / bed, title, out_png)
            panels.append(mpimg.imread(str(out_png)))

        fig, axes = plt.subplots(2, 2, figsize=(16, 10))
        for ax, img in zip(axes.flatten(), panels):
            ax.imshow(img)
            ax.axis("off")

        fig.suptitle("Metaplot of simulated region sets", fontsize=13)
        fig.tight_layout()
        OUT.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(OUT, dpi=150, bbox_inches="tight")
        print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
