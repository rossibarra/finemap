#!/usr/bin/env python3
"""Invert a UCSC chain file (target <-> query), like kent's chainSwap.

Used to turn the AnchorWave-derived v5 -> AGPv2 chain (data/v5v2.chain) into
an AGPv2 -> v5 chain (data/v2v5.chain) for CrossMap.

Chains with a '-' query strand have their query coordinates on the reverse
strand; after swapping, both sides are converted to reverse-strand coordinates
and the block order is reversed so the new target stays on '+'.
"""

import argparse
from pathlib import Path

ROOT = Path(__file__).parent.parent


def parse_args():
    parser = argparse.ArgumentParser(description="Swap target and query in a chain file.")
    parser.add_argument("chain_in", nargs="?", default=ROOT / "data/v5v2.chain", type=Path)
    parser.add_argument("chain_out", nargs="?", default=ROOT / "data/v2v5.chain", type=Path)
    return parser.parse_args()


def read_chains(path):
    header, lines = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("chain"):
                if header is not None:
                    yield header, lines
                header, lines = line.split(), []
            elif line:
                lines.append([int(x) for x in line.split()])
    if header is not None:
        yield header, lines


def swap(header, lines):
    _, score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, cid = header
    assert tStrand == "+", "target strand must be +"
    tSize, qSize = int(tSize), int(qSize)

    # Absolute (tStart, tEnd, qStart, qEnd) for each aligned block.
    blocks, t, q = [], int(tStart), int(qStart)
    for fields in lines:
        size = fields[0]
        blocks.append((t, t + size, q, q + size))
        if len(fields) == 3:
            t += size + fields[1]
            q += size + fields[2]
    assert blocks[-1][1] == int(tEnd) and blocks[-1][3] == int(qEnd)

    # Swap roles: new target = old query, new query = old target.
    blocks = [(qs, qe, ts, te) for ts, te, qs, qe in blocks]
    nt_name, nt_size, nq_name, nq_size = qName, qSize, tName, tSize
    if qStrand == "-":
        blocks = [
            (nt_size - te, nt_size - ts, nq_size - qe, nq_size - qs)
            for ts, te, qs, qe in reversed(blocks)
        ]

    out = [
        f"chain {score} {nt_name} {nt_size} + {blocks[0][0]} {blocks[-1][1]} "
        f"{nq_name} {nq_size} {qStrand} {blocks[0][2]} {blocks[-1][3]} {cid}"
    ]
    for (ts, te, qs, qe), nxt in zip(blocks, blocks[1:]):
        out.append(f"{te - ts}\t{nxt[0] - te}\t{nxt[2] - qe}")
    out.append(f"{blocks[-1][1] - blocks[-1][0]}")
    return "\n".join(out) + "\n\n"


def main():
    args = parse_args()
    n = 0
    with open(args.chain_out, "w") as out:
        for header, lines in read_chains(args.chain_in):
            out.write(swap(header, lines))
            n += 1
    print(f"Wrote {n} chains to {args.chain_out}")


if __name__ == "__main__":
    main()
