#!/usr/bin/env python3

from __future__ import annotations

import random
from collections import defaultdict
from pathlib import Path


JRI_PATH = Path("jri_v5.bed")
GFF_PATH = Path("v5.gff3")
FAI_PATH = Path("data/v5.fa.gz.fai")
README_PATH = Path("simulation_readme.md")

N_REGIONS = 100_000
FLANK_BP = 5_000
RANDOM_SEED = 25


def normalize_chrom(chrom: str) -> str:
    if chrom.startswith("Chr"):
        return "chr" + chrom[3:]
    if chrom.isdigit():
        return f"chr{chrom}"
    return chrom


def load_lengths() -> dict[str, int]:
    lengths: dict[str, int] = {}
    with FAI_PATH.open() as handle:
        for line in handle:
            seqid, length, *_ = line.rstrip("\n").split("\t")
            if seqid.startswith("chr") and seqid[3:].isdigit() and 1 <= int(seqid[3:]) <= 10:
                lengths[seqid] = int(length)
    return lengths


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    out = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start > out[-1][1]:
            out.append([start, end])
        elif end > out[-1][1]:
            out[-1][1] = end
    return [(start, end) for start, end in out]


def subtract_intervals(
    base: list[tuple[int, int]],
    subtract: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    if not base:
        return []
    if not subtract:
        return base[:]

    out: list[tuple[int, int]] = []
    j = 0
    for start, end in base:
        cur = start
        while j < len(subtract) and subtract[j][1] <= start:
            j += 1
        k = j
        while k < len(subtract) and subtract[k][0] < end:
            sub_start, sub_end = subtract[k]
            if sub_start > cur:
                out.append((cur, min(sub_start, end)))
            cur = max(cur, sub_end)
            if cur >= end:
                break
            k += 1
        if cur < end:
            out.append((cur, end))
    return [(s, e) for s, e in out if e > s]


def intersect_with_range(intervals: list[tuple[int, int]], lo: int, hi: int) -> list[tuple[int, int]]:
    out = []
    for start, end in intervals:
        s = max(start, lo)
        e = min(end, hi)
        if e > s:
            out.append((s, e))
    return out


def weighted_random_point(
    rng: random.Random,
    intervals: list[tuple[int, int]],
) -> int:
    total = sum(end - start for start, end in intervals)
    if total <= 0:
        raise ValueError("No available space to sample a point")
    offset = rng.randrange(total)
    for start, end in intervals:
        span = end - start
        if offset < span:
            return start + offset
        offset -= span
    raise RuntimeError("Failed to sample point")


def load_gene_partitions(lengths: dict[str, int]) -> tuple[
    dict[str, list[tuple[int, int]]],
    dict[str, list[tuple[int, int]]],
    dict[str, list[tuple[int, int]]],
]:
    five_prime: dict[str, list[tuple[int, int]]] = defaultdict(list)
    three_prime: dict[str, list[tuple[int, int]]] = defaultdict(list)
    body_core: dict[str, list[tuple[int, int]]] = defaultdict(list)

    with GFF_PATH.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = normalize_chrom(parts[0])
            if chrom not in lengths:
                continue
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            chrom_len = lengths[chrom]

            if strand == "+":
                five_prime[chrom].append((max(0, start - FLANK_BP), start))
                three_prime[chrom].append((end, min(chrom_len, end + FLANK_BP)))
                core_start = start + FLANK_BP
                core_end = end - FLANK_BP
            else:
                five_prime[chrom].append((end, min(chrom_len, end + FLANK_BP)))
                three_prime[chrom].append((max(0, start - FLANK_BP), start))
                core_start = start + FLANK_BP
                core_end = end - FLANK_BP

            if core_end > core_start:
                body_core[chrom].append((core_start, core_end))

    five_prime = {chrom: merge_intervals(vals) for chrom, vals in five_prime.items()}
    three_prime = {chrom: merge_intervals(vals) for chrom, vals in three_prime.items()}
    body_core = {chrom: merge_intervals(vals) for chrom, vals in body_core.items()}

    exclusive_five = {}
    exclusive_three = {}
    exclusive_body = {}
    for chrom in lengths:
        five = five_prime.get(chrom, [])
        three = three_prime.get(chrom, [])
        body = body_core.get(chrom, [])
        exclusive_five[chrom] = subtract_intervals(subtract_intervals(five, three), body)
        exclusive_three[chrom] = subtract_intervals(subtract_intervals(three, five), body)
        exclusive_body[chrom] = subtract_intervals(subtract_intervals(body, five), three)

    return exclusive_five, exclusive_three, exclusive_body


def load_length_distribution() -> list[int]:
    lengths = []
    with JRI_PATH.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            start = int(parts[1])
            end = int(parts[2])
            if end > start:
                lengths.append(end - start)
    return lengths


def sample_uniform_region(
    rng: random.Random,
    lengths_map: dict[str, int],
    length: int,
) -> tuple[str, int, int]:
    choices = []
    total = 0
    for chrom, chrom_len in lengths_map.items():
        available = chrom_len - length + 1
        if available > 0:
            total += available
            choices.append((chrom, total))
    pick = rng.randrange(total)
    for chrom, cumulative in choices:
        if pick < cumulative:
            chrom_len = lengths_map[chrom]
            start = rng.randrange(chrom_len - length + 1)
            return chrom, start, start + length
    raise RuntimeError("Uniform sampling failed")


def sample_point_from_annotation(
    rng: random.Random,
    intervals_by_chr: dict[str, list[tuple[int, int]]],
) -> tuple[str, int]:
    choices = []
    total = 0
    for chrom, intervals in intervals_by_chr.items():
        span = sum(end - start for start, end in intervals)
        if span > 0:
            total += span
            choices.append((chrom, total))
    pick = rng.randrange(total)
    for chrom, cumulative in choices:
        if pick < cumulative:
            point = weighted_random_point(rng, intervals_by_chr[chrom])
            return chrom, point
    raise RuntimeError("Point sampling failed")


def sample_enriched_region(
    rng: random.Random,
    intervals_by_chr: dict[str, list[tuple[int, int]]],
    lengths_map: dict[str, int],
    length: int,
) -> tuple[str, int, int]:
    chrom, point = sample_point_from_annotation(rng, intervals_by_chr)
    chrom_len = lengths_map[chrom]
    offset = rng.randrange(length)
    start = point - offset
    end = start + length
    if start < 0:
        start = 0
        end = length
    if end > chrom_len:
        end = chrom_len
        start = max(0, chrom_len - length)
    return chrom, start, end


def write_bed(path: Path, records: list[tuple[str, int, int, str]]) -> None:
    with path.open("w") as handle:
        for chrom, start, end, ident in records:
            handle.write(f"{chrom}\t{start}\t{end}\t{ident}\n")


def write_readme() -> None:
    README_PATH.write_text(
        "\n".join(
            [
                "# Simulated Example Region Sets",
                "",
                "These BED files each contain 100,000 simulated regions with lengths sampled from the empirical length distribution of `jri_v5.bed`.",
                "",
                "All coordinates are on B73 v5 chromosomes `chr1`-`chr10`.",
                "",
                "Sampling rules:",
                "",
                "- `example1.bed`: uniform along the genome. Region starts are sampled uniformly across all valid genomic positions for the chosen length.",
                "- `example2.bed`: enriched in 5' gene regions. One anchor point is sampled from strand-aware 5' 5 kb flanks, excluding positions that also fall in 3' flanks or the internal gene-body set below.",
                "- `example3.bed`: enriched in gene bodies. One anchor point is sampled from internal gene-body segments after removing the first and last 5 kb of each gene; this excludes both 5' and 3' 5 kb flanks.",
                "- `example4.bed`: enriched in 3' gene regions. One anchor point is sampled from strand-aware 3' 5 kb flanks, excluding positions that also fall in 5' flanks or the internal gene-body set.",
                "",
                "Notes:",
                "",
                "- Enrichment is defined by first drawing one point from the target annotation class, then assigning that point a random uniform position within the simulated region.",
                "- Gene annotations come from `v5.gff3`.",
                "- Chromosome lengths come from `data/v5.fa.gz.fai`.",
            ]
        )
        + "\n"
    )


def main() -> None:
    rng = random.Random(RANDOM_SEED)
    lengths_map = load_lengths()
    empirical_lengths = load_length_distribution()
    five_prime, three_prime, body_core = load_gene_partitions(lengths_map)

    outputs = {
        "example1": [],
        "example2": [],
        "example3": [],
        "example4": [],
    }

    for idx in range(1, N_REGIONS + 1):
        length = rng.choice(empirical_lengths)

        chrom, start, end = sample_uniform_region(rng, lengths_map, length)
        outputs["example1"].append((chrom, start, end, f"example1_{idx:06d}"))

        chrom, start, end = sample_enriched_region(rng, five_prime, lengths_map, length)
        outputs["example2"].append((chrom, start, end, f"example2_{idx:06d}"))

        chrom, start, end = sample_enriched_region(rng, body_core, lengths_map, length)
        outputs["example3"].append((chrom, start, end, f"example3_{idx:06d}"))

        chrom, start, end = sample_enriched_region(rng, three_prime, lengths_map, length)
        outputs["example4"].append((chrom, start, end, f"example4_{idx:06d}"))

    for name, records in outputs.items():
        write_bed(Path(f"{name}.bed"), records)
    write_readme()
    print("Wrote example1.bed through example4.bed and simulation_readme.md")


if __name__ == "__main__":
    main()
