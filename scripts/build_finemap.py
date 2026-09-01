#!/usr/bin/env python3
"""Derive finemap_v5.bed and per-chromosome HapMap files from jri_v5.bed."""

import sys
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).parent.parent
JRI = ROOT / "data/jri_v5.bed"
OGUT = ROOT / "data/ogut_fifthcM_map_agpv2.csv"
FAI = ROOT / "data/v5.fa.gz.fai"
FINEMAP_OUT = ROOT / "data/finemap_v5.bed"
HAPMAP_DIR = ROOT / "data/hapmap"


def sweep_line_weights(sub):
    """Return merged (start, end, weight) segments for one chromosome."""
    starts = pd.DataFrame({"pos": sub["start"].values, "delta": sub["weight"].values})
    ends = pd.DataFrame({"pos": sub["end"].values, "delta": -sub["weight"].values})
    events = pd.concat([starts, ends], ignore_index=True)
    pos_delta = events.groupby("pos")["delta"].sum().sort_index()

    positions = pos_delta.index.values
    cum_w = np.cumsum(pos_delta.values)

    segs = []
    for i in range(len(positions) - 1):
        w = cum_w[i]
        if w > 1e-15:
            segs.append([positions[i], positions[i + 1], w])

    if not segs:
        return []

    merged = [segs[0]]
    for s, e, w in segs[1:]:
        if abs(merged[-1][2] - w) < 1e-15:
            merged[-1][1] = e
        else:
            merged.append([s, e, w])
    return merged


def build_finemap(jri, chr_cM):
    records = []
    chroms = sorted(jri["chr"].unique(), key=lambda c: int(c.replace("Chr", "")))
    for chrom in chroms:
        sub = jri[jri["chr"] == chrom]
        merged = sweep_line_weights(sub)
        if not merged:
            continue

        target_cM = chr_cM.get(chrom)
        if target_cM is None:
            print(f"Warning: no Ogut cM for {chrom}", file=sys.stderr)
            continue

        total_mass = sum((e - s) * w for s, e, w in merged)
        scale = target_cM / total_mass

        cM_pos = 0.0
        for seg_s, seg_e, w in merged:
            cM_per_bp = w * scale
            cM_end = cM_pos + (seg_e - seg_s) * cM_per_bp

            if cM_end < cM_pos:
                raise FloatingPointError(
                    f"Decreasing cM on {chrom} at {seg_s}-{seg_e}: "
                    f"cM_start={cM_pos:.17g}, cM_end={cM_end:.17g}, "
                    f"weight={w:.17g}, scale={scale:.17g}"
                )
            cM_per_Mb = 0.0 if cM_end == cM_pos else cM_per_bp * 1e6

            records.append((chrom, seg_s, seg_e, cM_pos, cM_end, cM_per_Mb))
            cM_pos = cM_end
    return records


def build_hapmap(records, chr_lengths):
    hapmap_dir = HAPMAP_DIR
    hapmap_dir.mkdir(parents=True, exist_ok=True)

    chrom_records = {}
    for r in records:
        chrom_records.setdefault(r[0], []).append(r)

    for chrom, recs in chrom_records.items():
        num = chrom.replace("Chr", "")
        outfile = hapmap_dir / f"chr{num}.hapmap.tsv"

        breakpoints = {}  # pos -> cM

        def add_bp(pos, cM):
            if pos in breakpoints:
                if abs(breakpoints[pos] - cM) > 1e-9:
                    print(
                        f"Warning: conflicting cM at {chrom}:{pos}: "
                        f"{breakpoints[pos]} vs {cM}",
                        file=sys.stderr,
                    )
            else:
                breakpoints[pos] = cM

        first_start = recs[0][1]
        if first_start > 0:
            add_bp(0, 0.0)

        for _, seg_s, seg_e, cM_s, cM_e, _ in recs:
            add_bp(seg_s, cM_s)
            add_bp(seg_e, cM_e)

        chr_len = chr_lengths.get(chrom)
        if chr_len and chr_len not in breakpoints:
            last_cM = recs[-1][4]
            add_bp(chr_len, last_cM)

        sorted_pos = sorted(breakpoints)
        rows = []
        for i, pos in enumerate(sorted_pos):
            cM = breakpoints[pos]
            if i < len(sorted_pos) - 1:
                next_pos = sorted_pos[i + 1]
                next_cM = breakpoints[next_pos]
                rate = (next_cM - cM) / (next_pos - pos) * 1e6 if next_pos > pos else 0.0
            else:
                rate = 0.0
            rows.append((f"chr{num}", pos, rate, cM))

        with open(outfile, "w") as f:
            f.write("Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n")
            for chrom_name, pos, rate, cM in rows:
                f.write(f"{chrom_name}\t{pos}\t{rate:.10g}\t{cM:.10g}\n")

    print(f"Wrote hapmap files to {hapmap_dir}/")


def main():
    jri = pd.read_csv(JRI, sep="\t", header=None, names=["chr", "start", "end", "id", "src"])
    jri = jri[jri["end"] > jri["start"]].copy()
    jri["weight"] = 1.0 / (jri["end"] - jri["start"])

    ogut = pd.read_csv(OGUT)
    ogut["chr"] = "Chr" + ogut["chromosome"].astype(str)
    chr_cM = {
        chrom: grp["cM"].max() - grp["cM"].min()
        for chrom, grp in ogut.groupby("chr")
    }

    fai = pd.read_csv(FAI, sep="\t", header=None, usecols=[0, 1], names=["seq", "length"])
    fai["chr"] = fai["seq"].str.replace(r"^chr", "Chr", regex=True)
    chr_lengths = dict(zip(fai["chr"], fai["length"]))

    records = build_finemap(jri, chr_cM)

    with open(FINEMAP_OUT, "w") as f:
        for r in records:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\n")
    print(f"Wrote {len(records)} segments to {FINEMAP_OUT}")

    build_hapmap(records, chr_lengths)


if __name__ == "__main__":
    main()
