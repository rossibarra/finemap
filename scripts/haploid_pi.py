#!/usr/bin/env python3
"""Windowed nucleotide diversity (pi) from haploid all-sites VCFs.

The combined.chr*.all_sites.vcf.gz files carry one haploid allele per sample
(GT is a single character), so pixy -- which assumes diploid genotypes -- cannot
be used directly.  This computes the same estimator pixy does, site by site:

    diffs_site = (n^2 - sum_a c_a^2) / 2      pairs_site = n (n - 1) / 2

summed over sites in a window and divided, where n is the number of non-missing
haplotypes at that site and c_a the count of allele a.  Accumulating numerator
and denominator separately is what keeps missing data from biasing pi, and it
makes invariant sites count toward the denominator only -- which is the whole
reason an all-sites VCF is needed.

Output columns match pixy's pi table so downstream code is interchangeable.
"""

import argparse
import subprocess
import sys

import numpy as np

GT_TAIL = None  # bytes of genotype fields per line; derived from sample count


def read_populations(path):
    """Return {population: [sample, ...]} preserving file order."""
    pops = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                sys.exit(f"populations file needs 2 tab-separated columns: {line!r}")
            pops.setdefault(fields[1], []).append(fields[0])
    return pops


def read_windows(path, chrom_prefix):
    """Return {chrom: (starts, ends)} from a BED, 0-based half-open."""
    windows = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, start, end = line.split("\t")[:3]
            chrom = chrom_prefix(chrom)
            windows.setdefault(chrom, []).append((int(start), int(end)))
    out = {}
    for chrom, rows in windows.items():
        rows.sort()
        starts = np.array([r[0] for r in rows], dtype=np.int64)
        ends = np.array([r[1] for r in rows], dtype=np.int64)
        if starts[0] != 0 or not np.array_equal(starts[1:], ends[:-1]):
            sys.exit(f"{chrom}: windows must tile the chromosome from 0 without gaps")
        out[chrom] = (starts, ends)
    return out


def open_vcf(path):
    """Stream a (possibly non-BGZF) gzipped VCF through system gzip."""
    return subprocess.Popen(
        ["gzip", "-dc", path], stdout=subprocess.PIPE, bufsize=1 << 22
    )


def parse_header(stream):
    """Consume header lines; return the sample list."""
    for raw in stream:
        if raw.startswith(b"##"):
            continue
        if raw.startswith(b"#CHROM"):
            return raw.rstrip(b"\n").decode().split("\t")[9:]
        sys.exit("reached data lines without finding a #CHROM header")
    sys.exit("empty VCF")


def process_vcf(path, pop_index, windows, window_size, chunk_bytes, progress):
    """Accumulate diffs/comparisons/sites per population per window."""
    proc = open_vcf(path)
    stream = proc.stdout
    samples = parse_header(stream)
    n_samples = len(samples)
    tail_bytes = 2 * n_samples - 1  # single-char GT fields joined by tabs

    # Resolve sample names to column offsets once.
    offsets = {}
    for pop, names in pop_index.items():
        missing = [s for s in names if s not in samples]
        if missing:
            sys.exit(f"{path}: samples not in VCF: {', '.join(missing)}")
        offsets[pop] = np.array([samples.index(s) for s in names], dtype=np.int64)

    acc = {}
    chrom_seen = None
    n_lines = 0

    while True:
        lines = stream.readlines(chunk_bytes)
        if not lines:
            break

        chrom = lines[0].split(b"\t", 1)[0].decode()
        if chrom_seen is None:
            chrom_seen = chrom
            if chrom not in windows:
                sys.exit(f"{path}: chromosome {chrom} absent from the windows BED")
            n_win = len(windows[chrom][0])
            for pop in pop_index:
                acc[pop] = {
                    "diffs": np.zeros(n_win),
                    "pairs": np.zeros(n_win),
                    "sites": np.zeros(n_win),
                }
        elif chrom != chrom_seen:
            sys.exit(f"{path}: expected one chromosome per file, saw {chrom_seen} and {chrom}")

        # Genotype fields occupy a fixed-width tail, so slice them as a block
        # rather than splitting every line into 38 fields.
        try:
            tail = b"".join([ln[-(tail_bytes + 1):-1] for ln in lines])
            codes = (
                np.frombuffer(tail, dtype=np.uint8)
                .reshape(len(lines), tail_bytes)[:, ::2]
                .astype(np.int16)
                - 48
            )
        except ValueError:
            sys.exit(f"{path}: genotype fields are not the expected fixed width")
        # '.' (ASCII 46) becomes -2; treat anything below zero as missing.
        called = codes >= 0

        pos = np.fromiter(
            (int(ln.split(b"\t", 2)[1]) for ln in lines),
            dtype=np.int64,
            count=len(lines),
        )
        win = (pos - 1) // window_size  # VCF POS is 1-based, BED starts are 0-based
        np.clip(win, 0, n_win - 1, out=win)

        max_allele = int(codes.max()) if codes.size else 0
        for pop, idx in offsets.items():
            sub = codes[:, idx]
            sub_called = called[:, idx]
            n = sub_called.sum(axis=1).astype(np.int64)
            sumsq = np.zeros(len(lines), dtype=np.int64)
            for allele in range(max_allele + 1):
                c = (sub == allele).sum(axis=1).astype(np.int64)
                sumsq += c * c
            diffs = (n * n - sumsq) // 2
            pairs = n * (n - 1) // 2
            informative = n >= 2

            a = acc[pop]
            a["diffs"] += np.bincount(win, weights=diffs, minlength=n_win)
            a["pairs"] += np.bincount(win, weights=pairs, minlength=n_win)
            a["sites"] += np.bincount(win[informative], minlength=n_win)

        n_lines += len(lines)
        if progress:
            print(f"  {chrom}: {n_lines:,} sites", end="\r", file=sys.stderr)

    stream.close()
    if proc.wait() != 0:
        sys.exit(f"gzip failed on {path}")
    if progress:
        print(f"  {chrom_seen}: {n_lines:,} sites", file=sys.stderr)
    return chrom_seen, acc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", nargs="+", required=True, help="one gzipped VCF per chromosome")
    ap.add_argument("--populations", required=True, help="sample<TAB>population, no header")
    ap.add_argument("--windows", required=True, help="BED tiling each chromosome from 0")
    ap.add_argument("--window-size", type=int, default=100000)
    ap.add_argument("--out", required=True, help="output TSV (pixy pi format)")
    ap.add_argument("--lowercase-chrom", action="store_true",
                    help="map BED 'Chr1' to VCF 'chr1'")
    ap.add_argument("--chunk-bytes", type=int, default=1 << 24)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    prefix = (lambda c: c.replace("Chr", "chr")) if args.lowercase_chrom else (lambda c: c)
    pops = read_populations(args.populations)
    windows = read_windows(args.windows, prefix)

    if not args.quiet:
        sizes = ", ".join(f"{p} n={len(s)}" for p, s in pops.items())
        print(f"populations: {sizes}", file=sys.stderr)

    with open(args.out, "w") as out:
        out.write("pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi"
                  "\tno_sites\tcount_diffs\tcount_comparisons\n")
        for path in args.vcf:
            if not args.quiet:
                print(f"{path}", file=sys.stderr)
            chrom, acc = process_vcf(path, pops, windows, args.window_size,
                                     args.chunk_bytes, not args.quiet)
            starts, ends = windows[chrom]
            for pop in pops:
                a = acc[pop]
                with np.errstate(invalid="ignore", divide="ignore"):
                    pi = np.where(a["pairs"] > 0, a["diffs"] / a["pairs"], np.nan)
                for i in range(len(starts)):
                    value = "NA" if np.isnan(pi[i]) else f"{pi[i]:.8g}"
                    out.write(
                        f"{pop}\t{chrom}\t{starts[i] + 1}\t{ends[i]}\t{value}"
                        f"\t{int(a['sites'][i])}\t{int(a['diffs'][i])}"
                        f"\t{int(a['pairs'][i])}\n"
                    )

    if not args.quiet:
        print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
