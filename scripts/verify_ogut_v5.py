#!/usr/bin/env python3
"""Check lifted Ogut marker positions against the AGPv2 and v5 genome sequence.

For each marker in data/ogut_v5.csv, the AGPv2 sequence centred on the marker
(+-FLANK bp) is searched for within +-WINDOW bp of its v5 position on both
strands. A marker is supported when the best match has at most MAX_MISMATCH
mismatches and lies exactly at the stated v5 position.

Markers next to an indel between the assemblies fail that ungapped 101 bp
test even when correctly placed, so a marker is also supported when either
one-sided flank (the marker base plus FLANK bp to one side) matches at the
stated v5 position, on either strand, with at most MAX_FLANK_MISMATCH
mismatches. With --filter, unsupported markers are dropped and cM_norm is
recomputed.

Requires samtools and bgzipped, faidx-indexed genomes:
  data/B73_RefGen_v2.fa.gz  (https://download.maizegdb.org/B73_RefGen_v2/)
  data/Zm-B73-REFERENCE-NAM-5.0.fa.gz
"""

import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent
OGUT_V2 = ROOT / "data/ogut_fifthcM_map_agpv2.csv"
OGUT_V5 = ROOT / "data/ogut_v5.csv"
V2_FASTA = ROOT / "data/B73_RefGen_v2.fa.gz"
V5_FASTA = ROOT / "data/Zm-B73-REFERENCE-NAM-5.0.fa.gz"

FLANK = 50
WINDOW = 300
MAX_MISMATCH = 5
MAX_FLANK_MISMATCH = 2
RC = str.maketrans("ACGTN", "TGCAN")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Verify lifted Ogut v5 positions by sequence and optionally drop unsupported markers."
    )
    parser.add_argument("--filter", action="store_true",
                        help="Drop unsupported markers and overwrite data/ogut_v5.csv.")
    parser.add_argument("--report", type=Path,
                        help="Optional TSV with per-marker mismatches, offset and strand.")
    return parser.parse_args()


def fetch(fasta, regions, chunk=2000):
    seqs = []
    for i in range(0, len(regions), chunk):
        out = subprocess.run(["samtools", "faidx", str(fasta)] + regions[i:i + chunk],
                             capture_output=True, text=True, check=True).stdout
        for rec in out.split(">")[1:]:
            seqs.append("".join(rec.splitlines()[1:]).upper())
    return seqs


def best_match(query, target, centre):
    """Return (mismatches, offset of the query centre from target[centre], strand).

    The query has odd length with the marker at index FLANK, which is also the
    marker's index in its reverse complement, so one offset formula serves both strands.
    """
    length = len(query)
    if len(target) < length:
        return length, None, None
    windows = np.lib.stride_tricks.sliding_window_view(np.frombuffer(target.encode(), "S1"), length)
    result = (length + 1, None, None)
    for strand, seq in (("+", query), ("-", query.translate(RC)[::-1])):
        mismatches = (windows != np.frombuffer(seq.encode(), "S1")).sum(axis=1)
        i = int(mismatches.argmin())
        if mismatches[i] < result[0]:
            result = (int(mismatches[i]), i + FLANK - centre, strand)
    return result


def flank_mismatches(query, target, centre):
    """Fewest mismatches of a one-sided flank placed exactly at target[centre], either strand."""
    best = len(query)
    for seq in (query, query.translate(RC)[::-1]):
        for q, lo in ((seq[:FLANK + 1], centre - FLANK), (seq[FLANK:], centre)):
            t = target[lo:lo + FLANK + 1] if lo >= 0 else ""
            if len(t) == len(q):
                best = min(best, sum(a != b for a, b in zip(q, t)))
    return best


def main():
    args = parse_args()
    ogut = pd.read_csv(OGUT_V5)
    v2 = pd.read_csv(OGUT_V2)[["SNP_newID", "position"]]
    markers = ogut.merge(v2, on="SNP_newID", how="left")
    num = markers["chr"].str.replace("Chr", "")

    queries = fetch(V2_FASTA, [f"Chr{c}:{p - FLANK}-{p + FLANK}" for c, p in zip(num, markers["position"])])
    targets = fetch(V5_FASTA, [f"chr{c}:{max(1, p - WINDOW)}-{p + WINDOW}" for c, p in zip(num, markers["pos_v5"])])
    # Index of the v5 position in its window (less than WINDOW near a chromosome start)
    centres = [p - max(1, p - WINDOW) for p in markers["pos_v5"]]
    matches = [best_match(q, t, c) for q, t, c in zip(queries, targets, centres)]
    markers["mismatches"], markers["offset"], markers["strand"] = zip(*matches)
    markers["flank_mismatches"] = [flank_mismatches(q, t, c) for q, t, c in zip(queries, targets, centres)]
    full = (markers["mismatches"] <= MAX_MISMATCH) & (markers["offset"] == 0)
    flank = markers["flank_mismatches"] <= MAX_FLANK_MISMATCH
    markers["supported"] = full | flank

    n_bad = (~markers["supported"]).sum()
    print(f"{markers['supported'].sum()} of {len(markers)} markers supported by sequence "
          f"({full.sum()} full window, {(flank & ~full).sum()} one flank only); {n_bad} unsupported")
    if args.report:
        markers.to_csv(args.report, sep="\t", index=False)
        print(f"Wrote {args.report}")

    if args.filter:
        kept = markers.loc[markers["supported"], ["chr", "pos_v5", "SNP_newID", "cM"]].copy()
        kept["cM_norm"] = kept.groupby("chr")["cM"].transform(lambda x: x - x.min())
        kept.to_csv(OGUT_V5, index=False)
        print(f"Wrote {len(kept)} markers to {OGUT_V5}")


if __name__ == "__main__":
    main()
