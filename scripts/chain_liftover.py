"""Single-base UCSC chain mapping with alignment identity preserved.

Coordinates in this module are zero-based. Each hit retains its chain and strand;
ambiguous hits are deliberately not resolved by alignment score.
"""

from collections import defaultdict
from dataclasses import dataclass
import heapq

from swap_chain import read_chains


@dataclass(frozen=True)
class Hit:
    chrom: str
    position: int
    strand: str
    chain: str


def map_points(chain_path, points):
    """Map {source chromosome: iterable of positions} to all hits per point.

    A sweep over aligned blocks handles overlapping chains without choosing an
    arbitrary mapping. Memory scales with blocks and unique queried endpoints.
    """
    blocks = defaultdict(list)
    for header, lines in read_chains(chain_path):
        _, _, chrom, size, strand, start, end, dest, dsize, dstrand, ds, de, cid = header
        if strand != "+" or dstrand not in {"+", "-"}:
            raise ValueError(f"Unsupported strand in chain {cid}")
        start, end, size, dsize, ds, de = map(int, (start, end, size, dsize, ds, de))
        t, q = start, ds
        if not (0 <= start < end <= size and 0 <= ds < de <= dsize):
            raise ValueError(f"Invalid bounds in chain {cid}")
        for fields in lines:
            width = fields[0]
            if width < 0 or len(fields) not in (1, 3):
                raise ValueError(f"Invalid block in chain {cid}")
            # AnchorWave chains contain zero-size blocks between adjacent gaps.
            # Advance across their gaps but never treat them as aligned bases.
            if width:
                blocks[chrom].append((t, t + width, dest, q, dsize, dstrand, cid))
            t += width
            q += width
            if len(fields) == 3:
                if min(fields[1:]) < 0:
                    raise ValueError(f"Negative chain gap in {cid}")
                t += fields[1]
                q += fields[2]
        if (t, q) != (end, de):
            raise ValueError(f"Blocks do not span header in chain {cid}")

    mapped = {}
    for chrom, positions in points.items():
        ordered = sorted(blocks.get(chrom, []))
        active, endings, index = {}, [], 0
        for position in sorted(set(positions)):
            while index < len(ordered) and ordered[index][0] <= position:
                active[index] = ordered[index]
                heapq.heappush(endings, (ordered[index][1], index))
                index += 1
            while endings and endings[0][0] <= position:
                active.pop(heapq.heappop(endings)[1], None)
            hits = []
            for start, _, dest, ds, dsize, strand, cid in active.values():
                p = ds + position - start
                if strand == "-":
                    p = dsize - 1 - p
                hits.append(Hit(dest, p, strand, cid))
            mapped[chrom, position] = hits
    return mapped
