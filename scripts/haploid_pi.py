#!/usr/bin/env python3
"""Windowed nucleotide diversity (pi) from haploid all-sites VCFs.

The combined.chr*.all_sites.vcf.gz files carry one haploid allele per sample,
so pixy -- which assumes diploid genotypes -- cannot
be used directly.  This computes the same estimator pixy does, site by site:

    diffs_site = (n^2 - sum_a c_a^2) / 2      pairs_site = n (n - 1) / 2

summed over sites in a window and divided, where n is the number of non-missing
haplotypes at that site and c_a the count of allele a.  Accumulating numerator
and denominator separately is what keeps missing data from biasing pi, and it
makes invariant sites count toward the denominator only -- which is the whole
reason an all-sites VCF is needed.

Output columns match pixy's pi table so downstream code is interchangeable.

Supported input: FORMAT must contain exactly one GT key (GT, GT:DP, DP:GT, ...)
and every GT must be a single haploid allele index or '.'.  Records whose
FORMAT is exactly GT with one-character sample fields take a vectorised fast
path; anything else is parsed field by field.  Diploid/polyploid calls
(0/0, 0|1, ./.), alleles not listed in ALT, and records without GT are
rejected with an error.  Run ``haploid_pi.py --self-test`` to check parsing.
"""

import argparse
import subprocess
import sys

import numpy as np

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
        if (starts[0] != 0 or np.any(ends <= starts)
                or not np.array_equal(starts[1:], ends[:-1])):
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


def load_restrict(path, key):
    """Return {chrom: sorted 1-based position array} from an .npz of '<chrom>_<key>' arrays."""
    if path is None:
        return None
    data = np.load(path)
    suffix = f"_{key}"
    out = {}
    for name in data.files:
        if name.endswith(suffix):
            arr = np.asarray(data[name], dtype=np.int64)
            out[name[: -len(suffix)]] = np.sort(arr)
    if not out:
        sys.exit(f"{path}: no arrays ending in {suffix!r}")
    return out


def parse_records_fast(lines, n_samples, chrom_seen, path):
    """Vectorised parser for the common case; return None to fall back.

    Applies only when every record in the chunk has FORMAT exactly ``GT`` and
    every sample field is one character ('0'-'9' or '.'), i.e. the layout of
    the combined.chr*.all_sites.vcf.gz files.  Anything else -- other FORMAT
    keys, multi-character or diploid GTs, CRLF, a missing final newline --
    returns None so parse_records() handles (or rejects) it field by field.
    """
    buf = b"".join(lines)
    if b"\r" in buf:
        return None
    arr = np.frombuffer(buf, dtype=np.uint8)
    nl = np.flatnonzero(arr == 10)
    if len(nl) != len(lines) or nl[-1] != len(arr) - 1:
        return None
    # Exactly 8 + n_samples tabs on every line: the right total, with each
    # line's share of the (sorted) tab positions lying between its newlines.
    tabs = np.flatnonzero(arr == 9)
    if tabs.size != len(nl) * (8 + n_samples):
        return None
    tabs = tabs.reshape(len(nl), 8 + n_samples)
    if np.any(tabs[:, -1] > nl) or np.any(tabs[1:, 0] < nl[:-1]):
        return None
    # Each line must end "\tGT\t" + n_samples one-byte fields joined by tabs.
    width = 2 * n_samples + 3
    if np.any(np.diff(nl, prepend=-1) <= width):
        return None
    block = arr[nl[:, None] + np.arange(-width, 0)]
    if not np.all(block[:, :4] == np.frombuffer(b"\tGT\t", dtype=np.uint8)):
        return None
    gt = block[:, 4::2]
    if n_samples > 1 and not np.all(block[:, 5::2] == 9):
        return None
    digit = (gt >= 48) & (gt <= 57)
    if not np.all(digit | (gt == 46)):
        return None
    codes = np.where(digit, gt.astype(np.int64) - 48, -1)

    parts = [ln.split(b"\t", 2) for ln in lines]
    chroms = {p[0] for p in parts}
    first = parts[0][0].decode()
    if chrom_seen is None:
        chrom_seen = first
    if len(chroms) != 1 or first != chrom_seen:
        return None  # let parse_records name the offending chromosome
    try:
        positions = np.fromiter((int(p[1]) for p in parts), dtype=np.int64,
                                count=len(parts))
    except ValueError:
        return None
    if np.any(positions < 1):
        return None
    # Allele indices must exist in ALT; only rows carrying a non-ref call matter.
    row_max = codes.max(axis=1)
    rows = np.flatnonzero(row_max > 0)
    if rows.size:
        alts = [lines[r].split(b"\t", 5)[4] for r in rows]
        n_alt = np.array([0 if a == b"." else a.count(b",") + 1 for a in alts])
        if np.any(row_max[rows] > n_alt):
            return None
    return chrom_seen, positions, codes


def parse_records(lines, n_samples, chrom_seen, path):
    """Parse haploid GT by FORMAT; reject malformed or mixed-chromosome records.

    Supported: FORMAT containing exactly one GT key (in any position, e.g. GT,
    GT:DP, DP:GT); per-sample GT values that are a single haploid allele index
    (0, 1, ..., up to the number of ALT alleles) or '.'; a bare '.' sample
    field.  Diploid or polyploid calls (0/0, 0|1, ./.) and alleles not
    present in ALT are rejected with an error rather than guessed at.
    """
    fast = parse_records_fast(lines, n_samples, chrom_seen, path)
    if fast is not None:
        return fast
    positions = np.empty(len(lines), dtype=np.int64)
    codes = np.full((len(lines), n_samples), -1, dtype=np.int64)
    for row, raw in enumerate(lines):
        fields = raw.rstrip(b"\r\n").split(b"\t")
        if len(fields) != 9 + n_samples:
            sys.exit(f"{path}: expected {9 + n_samples} VCF columns")
        chrom = fields[0].decode()
        if chrom_seen is None:
            chrom_seen = chrom
        elif chrom != chrom_seen:
            sys.exit(f"{path}: expected one chromosome per file, saw {chrom_seen} and {chrom}")
        try:
            pos = int(fields[1])
        except ValueError:
            sys.exit(f"{path}: invalid VCF position {fields[1]!r}")
        if pos < 1:
            sys.exit(f"{path}: VCF position must be positive: {pos}")
        positions[row] = pos
        fmt = fields[8].split(b":")
        if fmt.count(b"GT") != 1:
            sys.exit(f"{path}: {chrom}:{pos}: FORMAT must contain exactly one GT")
        gt_index = fmt.index(b"GT")
        n_alt = 0 if fields[4] == b"." else len(fields[4].split(b","))
        for col, sample in enumerate(fields[9:]):
            if sample == b".":
                continue
            values = sample.split(b":")
            if len(values) > len(fmt) or len(values) <= gt_index:
                sys.exit(f"{path}: {chrom}:{pos}: sample fields do not match FORMAT")
            gt = values[gt_index]
            if gt == b".":
                continue
            if b"/" in gt or b"|" in gt:
                sys.exit(f"{path}: {chrom}:{pos}: diploid/polyploid GT {gt!r} is not "
                         "supported; this estimator expects one haploid allele per sample")
            if not gt.isdigit() or int(gt) > n_alt:
                sys.exit(f"{path}: {chrom}:{pos}: expected haploid GT allele in 0..{n_alt}, got {gt!r}")
            codes[row, col] = int(gt)
    return chrom_seen, positions, codes


def process_vcf(path, pop_index, windows, window_size, chunk_bytes, progress,
                restrict=None):
    """Accumulate diffs/comparisons/sites per population per window."""
    proc = open_vcf(path)
    try:
        result = accumulate_stream(proc.stdout, path, pop_index, windows,
                                   chunk_bytes, progress, restrict)
        if proc.wait() != 0:
            sys.exit(f"gzip failed on {path}")
        return result
    finally:
        proc.stdout.close()
        if proc.poll() is None:
            proc.terminate()
        proc.wait()


def accumulate_stream(stream, path, pop_index, windows, chunk_bytes, progress,
                      restrict=None):
    samples = parse_header(stream)
    n_samples = len(samples)
    if not samples or len(set(samples)) != n_samples:
        sys.exit(f"{path}: VCF sample names must be nonempty and unique")

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

        chrom, pos, codes = parse_records(lines, n_samples, chrom_seen, path)
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
        called = codes >= 0
        starts, ends = windows[chrom]
        win = np.searchsorted(starts, pos - 1, side="right") - 1
        if np.any(win < 0) or np.any(pos > ends[win]):
            sys.exit(f"{path}: VCF position outside windows BED")
        if restrict is not None:
            # Membership test against the sorted position list for this chromosome,
            # rather than a per-chromosome boolean mask of ~300M elements.
            target = restrict.get(chrom_seen)
            if target is None or target.size == 0:
                keep = np.zeros(len(pos), dtype=bool)
            else:
                j = np.searchsorted(target, pos)
                np.clip(j, 0, target.size - 1, out=j)
                keep = target[j] == pos
            if not keep.any():
                n_lines += len(lines)
                continue
            pos = pos[keep]
            codes = codes[keep]
            called = called[keep]
            win = win[keep]

        max_allele = int(codes.max()) if codes.size else 0
        for pop, idx in offsets.items():
            sub = codes[:, idx]
            sub_called = called[:, idx]
            n = sub_called.sum(axis=1).astype(np.int64)
            sumsq = np.zeros(len(pos), dtype=np.int64)
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

    if progress:
        print(f"  {chrom_seen}: {n_lines:,} sites", file=sys.stderr)
    if chrom_seen is None:
        sys.exit(f"{path}: VCF contains no records")
    return chrom_seen, acc


def self_test():
    """Check GT parsing on synthetic records; exits non-zero on failure."""
    import io

    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n"
    windows = {"chr1": (np.array([0]), np.array([100]))}

    def run(records):
        stream = io.BytesIO((header + records).encode())
        acc = accumulate_stream(stream, "self-test", {"p": ["a", "b", "c"]},
                                windows, 1 << 20, False)[1]["p"]
        return int(acc["diffs"][0]), int(acc["pairs"][0])

    def rec(pos, fmt, *samples, alt="C"):
        return f"chr1\t{pos}\t.\tA\t{alt}\t.\tPASS\t.\t{fmt}\t" + "\t".join(samples) + "\n"

    checks = [
        # GT:DP must not read DP digits as alleles (the old fixed-width bug).
        ("GT:DP invariant", rec(1, "GT:DP", "0:8", "0:9", "0:7"), (0, 3)),
        ("GT:DP variant", rec(1, "GT:DP", "0:8", "1:9", "0:7"), (2, 3)),
        ("GT not first", rec(1, "DP:GT", "8:0", "9:1", "7:1"), (2, 3)),
        # Haploid missing calls, as '.' GT or a bare '.' sample, drop out of n.
        ("haploid missing", rec(1, "GT", "0", ".", "1"), (1, 1)),
        ("missing with DP", rec(1, "GT:DP", ".:3", "1:4", "."), (0, 0)),
        ("fast path", rec(1, "GT", "0", "1", "1") + rec(2, "GT", "0", "0", "."), (2, 4)),
        ("mixed chunk", rec(1, "GT", "0", "1", "1") + rec(2, "GT:DP", "0:1", "0:2", "1:3"),
         (4, 6)),
    ]
    failed = 0
    for name, records, want in checks:
        got = run(records)
        if got != want:
            print(f"FAIL {name}: (diffs, pairs) {got} != {want}", file=sys.stderr)
            failed += 1
    rejected = [
        ("diploid GT", rec(1, "GT", "0/0", "0|1", "./.")),
        ("diploid GT:DP", rec(1, "GT:DP", "0/1:5", "0:3", "0:3")),
        ("no GT key", rec(1, "DP", "8", "9", "7")),
        ("allele absent from ALT", rec(1, "GT", "0", "2", "0")),
        ("invariant site with ALT call", rec(1, "GT", "0", "1", "0", alt=".")),
    ]
    for name, records in rejected:
        try:
            run(records)
        except SystemExit:
            continue
        print(f"FAIL {name}: record was accepted", file=sys.stderr)
        failed += 1
    if failed:
        sys.exit(f"haploid_pi self-test: {failed} check(s) failed")
    print("haploid_pi self-test passed")


def main():
    if "--self-test" in sys.argv[1:]:
        self_test()
        return
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", nargs="+", required=True, help="one gzipped VCF per chromosome")
    ap.add_argument("--populations", required=True, help="sample<TAB>population, no header")
    ap.add_argument("--windows", required=True, help="BED tiling each chromosome from 0")
    ap.add_argument("--window-size", type=int, default=100000,
                    help="legacy compatibility option; windows are defined by --windows")
    ap.add_argument("--out", required=True, help="output TSV (pixy pi format)")
    ap.add_argument("--lowercase-chrom", action="store_true",
                    help="map BED 'Chr1' to VCF 'chr1'")
    ap.add_argument("--restrict", help="npz of '<chrom>_<key>' 1-based position arrays; "
                                       "only those sites contribute to pi")
    ap.add_argument("--restrict-key", default="4D",
                    help="array suffix inside --restrict (default 4D)")
    ap.add_argument("--chunk-bytes", type=int, default=1 << 24)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    prefix = (lambda c: c.replace("Chr", "chr")) if args.lowercase_chrom else (lambda c: c)
    pops = read_populations(args.populations)
    windows = read_windows(args.windows, prefix)
    restrict = load_restrict(args.restrict, args.restrict_key)
    if restrict is not None and not args.quiet:
        total = sum(a.size for a in restrict.values())
        print(f"restricting to {total:,} {args.restrict_key} sites", file=sys.stderr)

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
                                     args.chunk_bytes, not args.quiet, restrict)
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
