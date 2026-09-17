#!/usr/bin/env python3
"""Annotate 0-fold and 4-fold degenerate sites from a reference FASTA + GFF3.

A CDS position is 4-fold degenerate when all three alternative bases leave the
encoded amino acid unchanged, and 0-fold degenerate when all three change it.
Everything else (2-fold, 3-fold) is discarded, which is what makes 4D a clean
synonymous proxy and 0D a clean nonsynonymous proxy: mutations at those sites
are unambiguously one or the other, without needing to know the derived allele
or worry about which of several possible changes occurred.

One transcript per gene is used -- the one flagged `canonical_transcript=1` in
the GFF3, else the one with the longest CDS -- so that alternative transcripts
of the same gene do not contribute a position twice.  Positions claimed by more
than one gene with conflicting classes (overlapping genes on opposite strands)
are dropped.

Writes a compact .npz holding, per chromosome, sorted int32 arrays of 1-based
positions: `<chrom>_0D` and `<chrom>_4D`.
"""

import argparse
import gzip
import subprocess
import sys
from collections import defaultdict

import numpy as np

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
BASES = "ACGT"  # code 0..3; anything else (N, IUPAC) becomes 4


def build_degeneracy_table():
    """(64, 3) int8: 4 for 4-fold, 0 for 0-fold, -1 otherwise; stops are -1."""
    table = np.full((64, 3), -1, dtype=np.int8)
    for i in range(64):
        codon = BASES[i // 16] + BASES[(i // 4) % 4] + BASES[i % 4]
        aa = CODON_TABLE[codon]
        if aa == "*":
            continue  # stop codons are excluded outright
        for p in range(3):
            same = 0
            for b in BASES:
                if b == codon[p]:
                    continue
                alt = codon[:p] + b + codon[p + 1:]
                if CODON_TABLE[alt] == aa:
                    same += 1
            if same == 3:
                table[i, p] = 4
            elif same == 0:
                table[i, p] = 0
    return table


def base_codes():
    """256-entry uint8 lookup: A/C/G/T (either case) -> 0..3, everything else 4."""
    lut = np.full(256, 4, dtype=np.uint8)
    for code, b in enumerate(BASES):
        lut[ord(b)] = code
        lut[ord(b.lower())] = code
    return lut


def read_gff(path, chroms):
    """Return {gene: {transcript: (chrom, strand, [(start, end, phase), ...], canonical)}}."""
    opener = gzip.open if str(path).endswith(".gz") else open
    tx_meta = {}      # transcript -> (gene, chrom, strand, canonical)
    tx_cds = defaultdict(list)
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[0] not in chroms:
                continue
            if f[2] == "mRNA":
                attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                if attrs.get("biotype") != "protein_coding":
                    continue
                tx_meta[attrs["ID"]] = (attrs.get("Parent", attrs["ID"]), f[0], f[6],
                                        attrs.get("canonical_transcript") == "1")
            elif f[2] == "CDS":
                attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                parent = attrs.get("Parent")
                if parent is None:
                    continue
                phase = 0 if f[7] == "." else int(f[7])
                tx_cds[parent].append((int(f[3]), int(f[4]), phase))
    return tx_meta, tx_cds


def pick_transcripts(tx_meta, tx_cds):
    """One transcript per gene: canonical if flagged, else longest CDS."""
    best = {}
    for tx, (gene, chrom, strand, canonical) in tx_meta.items():
        cds = tx_cds.get(tx)
        if not cds:
            continue
        length = sum(e - s + 1 for s, e, _ in cds)
        key = (1 if canonical else 0, length)
        if gene not in best or key > best[gene][0]:
            best[gene] = (key, tx, chrom, strand, cds)
    return [(tx, chrom, strand, cds) for _, (_, tx, chrom, strand, cds) in best.items()]


def iter_fasta(path, wanted):
    """Yield (name, uint8 array of ASCII bytes) for wanted sequences, streaming."""
    proc = subprocess.Popen(["gzip", "-dc", str(path)], stdout=subprocess.PIPE,
                            bufsize=1 << 22)
    name, chunks, keep = None, [], False
    for raw in proc.stdout:
        if raw.startswith(b">"):
            if keep and name is not None:
                yield name, np.frombuffer(b"".join(chunks), dtype=np.uint8)
            name = raw[1:].split()[0].decode()
            keep = name in wanted
            chunks = []
        elif keep:
            chunks.append(raw.rstrip(b"\n"))
    if keep and name is not None:
        yield name, np.frombuffer(b"".join(chunks), dtype=np.uint8)
    proc.stdout.close()
    proc.wait()


def classify_chromosome(seq_codes, transcripts, deg_table):
    """Return (pos0, pos4, codonpos0, codonpos4) as 1-based positions, unsorted."""
    out_pos = {0: [], 4: []}
    out_cp = {0: [], 4: []}
    n = len(seq_codes)
    for _tx, _chrom, strand, cds in transcripts:
        cds = sorted(cds)
        if cds[0][0] < 1 or cds[-1][1] > n:
            continue
        # Genomic positions in transcript (5'->3') order.
        pieces = [np.arange(s, e + 1, dtype=np.int64) for s, e, _ in cds]
        if strand == "-":
            pieces = [p[::-1] for p in pieces[::-1]]
            first_phase = cds[-1][2]
        else:
            first_phase = cds[0][2]
        pos = np.concatenate(pieces)
        bases = seq_codes[pos - 1]
        if strand == "-":
            bases = np.where(bases < 4, 3 - bases, 4).astype(np.uint8)
        # Drop the leading partial codon indicated by phase, then the trailing remainder.
        if first_phase:
            pos = pos[first_phase:]
            bases = bases[first_phase:]
        ncod = len(pos) // 3
        if ncod == 0:
            continue
        pos = pos[: ncod * 3].reshape(ncod, 3)
        bases = bases[: ncod * 3].reshape(ncod, 3).astype(np.int64)
        good = (bases < 4).all(axis=1)
        if not good.any():
            continue
        pos, bases = pos[good], bases[good]
        idx = bases[:, 0] * 16 + bases[:, 1] * 4 + bases[:, 2]
        cls = deg_table[idx]  # (ncod, 3)
        for fold in (0, 4):
            rows, cols = np.nonzero(cls == fold)
            out_pos[fold].append(pos[rows, cols])
            out_cp[fold].append(cols.astype(np.int8) + 1)
    def cat(lst, dtype):
        return np.concatenate(lst).astype(dtype) if lst else np.empty(0, dtype=dtype)
    return (cat(out_pos[0], np.int64), cat(out_pos[4], np.int64),
            cat(out_cp[0], np.int8), cat(out_cp[4], np.int8))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="reference FASTA (gzipped ok)")
    ap.add_argument("--gff", required=True, help="full GFF3 with CDS features (gzipped ok)")
    ap.add_argument("--chroms", nargs="+",
                    default=[f"chr{i}" for i in range(1, 11)])
    ap.add_argument("--out", required=True, help="output .npz")
    args = ap.parse_args()

    chroms = list(args.chroms)
    deg_table = build_degeneracy_table()
    lut = base_codes()

    print(f"reading {args.gff}", file=sys.stderr)
    tx_meta, tx_cds = read_gff(args.gff, set(chroms))
    chosen = pick_transcripts(tx_meta, tx_cds)
    n_canon = sum(1 for tx, *_ in chosen if tx_meta[tx][3])
    print(f"  {len(tx_meta):,} protein-coding transcripts -> {len(chosen):,} genes "
          f"({n_canon:,} canonical-flagged)", file=sys.stderr)

    by_chrom = defaultdict(list)
    for rec in chosen:
        by_chrom[rec[1]].append(rec)

    arrays = {}
    totals = {0: 0, 4: 0}
    codon_pos_counts = {0: np.zeros(4, dtype=np.int64), 4: np.zeros(4, dtype=np.int64)}
    dropped_conflict = 0

    for name, seq in iter_fasta(args.fasta, set(chroms)):
        codes = lut[seq]
        p0, p4, c0, c4 = classify_chromosome(codes, by_chrom.get(name, []), deg_table)
        # Collapse duplicates from overlapping genes; drop positions claimed by both.
        u0, i0 = np.unique(p0, return_index=True)
        u4, i4 = np.unique(p4, return_index=True)
        conflict = np.intersect1d(u0, u4, assume_unique=True)
        if conflict.size:
            dropped_conflict += conflict.size
            keep0 = ~np.isin(u0, conflict)
            keep4 = ~np.isin(u4, conflict)
            u0, i0 = u0[keep0], i0[keep0]
            u4, i4 = u4[keep4], i4[keep4]
        arrays[f"{name}_0D"] = u0.astype(np.int32)
        arrays[f"{name}_4D"] = u4.astype(np.int32)
        totals[0] += len(u0)
        totals[4] += len(u4)
        codon_pos_counts[0] += np.bincount(c0[i0], minlength=4)
        codon_pos_counts[4] += np.bincount(c4[i4], minlength=4)
        print(f"  {name}: {len(u0):,} 0D, {len(u4):,} 4D", file=sys.stderr)

    np.savez_compressed(args.out, **arrays)
    print(f"\nwrote {args.out}", file=sys.stderr)
    print(f"genome-wide: {totals[0]:,} 0D sites, {totals[4]:,} 4D sites "
          f"({dropped_conflict:,} dropped for conflicting overlap)", file=sys.stderr)
    for fold in (0, 4):
        c = codon_pos_counts[fold][1:]
        frac = c / max(c.sum(), 1)
        print(f"  {fold}D codon-position split: "
              + ", ".join(f"pos{i+1} {frac[i]*100:.2f}%" for i in range(3)), file=sys.stderr)


if __name__ == "__main__":
    main()
