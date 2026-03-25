#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path


GFF_PATH = Path("v5.gff3")
FAI_PATH = Path("v5.fa.gz.fai")
UPSTREAM_OUTPUT = Path("results/v5_genes_5prime_upstream_2kb.bed")
DOWNSTREAM_OUTPUT = Path("results/v5_genes_3prime_downstream_2kb.bed")
FLANK_BP = 2000


def load_lengths() -> dict[str, int]:
    lengths: dict[str, int] = {}
    with FAI_PATH.open() as handle:
        for line in handle:
            seqid, length, *_rest = line.rstrip("\n").split("\t")
            lengths[seqid] = int(length)
    return lengths


def parse_gene_id(attributes: str) -> str:
    for field in attributes.split(";"):
        if field.startswith("ID="):
            return field[3:]
    return "unknown_gene"


def gff_interval_to_bed(start_1based: int, end_1based: int) -> tuple[int, int]:
    return start_1based - 1, end_1based


def main() -> None:
    seq_lengths = load_lengths()
    UPSTREAM_OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    with UPSTREAM_OUTPUT.open("w", newline="") as upstream_handle, DOWNSTREAM_OUTPUT.open("w", newline="") as downstream_handle:
        upstream_writer = csv.writer(upstream_handle, delimiter="\t")
        downstream_writer = csv.writer(downstream_handle, delimiter="\t")

        with GFF_PATH.open() as gff_handle:
            for line in gff_handle:
                if line.startswith("#"):
                    continue

                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue

                seqid, _source, _feature_type, start_str, end_str, _score, strand, _phase, attributes = parts
                if seqid not in seq_lengths:
                    continue

                chrom_length = seq_lengths[seqid]
                start_1based = int(start_str)
                end_1based = int(end_str)
                gene_id = parse_gene_id(attributes)

                gene_bed_start, gene_bed_end = gff_interval_to_bed(start_1based, end_1based)

                if strand == "+":
                    upstream_start = max(0, gene_bed_start - FLANK_BP)
                    upstream_end = gene_bed_start
                    downstream_start = gene_bed_end
                    downstream_end = min(chrom_length, gene_bed_end + FLANK_BP)
                elif strand == "-":
                    upstream_start = gene_bed_end
                    upstream_end = min(chrom_length, gene_bed_end + FLANK_BP)
                    downstream_start = max(0, gene_bed_start - FLANK_BP)
                    downstream_end = gene_bed_start
                else:
                    continue

                if upstream_start < upstream_end:
                    upstream_writer.writerow([seqid, upstream_start, upstream_end, gene_id, strand])

                if downstream_start < downstream_end:
                    downstream_writer.writerow([seqid, downstream_start, downstream_end, gene_id, strand])

    print(f"Wrote {UPSTREAM_OUTPUT} and {DOWNSTREAM_OUTPUT}")


if __name__ == "__main__":
    main()
