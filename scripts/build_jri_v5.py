#!/usr/bin/env python3
"""Lift the four crossover interval sources to B73 v5 and build data/jri_v5.bed.

Runs README Steps 2 and 3:

1. Convert each source to BED with the chromosome names its chain expects
   (AGPv2 chain: Chr1..Chr10; AGPv4 chain: 1..10).
2. Source coordinates are 1-based SNP positions, so an interval spans the
   1-based closed range [start, end], i.e. BED [start-1, end). Split it into
   two 1-bp endpoint markers, BED [start-1, start) and [end-1, end), and lift
   them with CrossMap.
3. Keep an interval only when both endpoints lift to the same chromosome and
   that chromosome is 1..10; the v5 interval is [min lifted start, max lifted end).
4. Normalise chromosome names to Chr1..Chr10, assign IDs and write
     data/jri_rm_euro_v5.bed                                   (5 cols, RMv2_/EUROv2_ IDs)
     data/xo_Combined_LR13_LR14_parents2_v5.bed                (4 cols)
     data/xo_ZeaGBSv27raw_RareAllelesC2TeoCurated_depth_v5.bed (4 cols)
     data/jri_v5.bed                                           (5 cols: chr, start, end, sample, id)
5. Validate chromosome names, schema, IDs, coordinates and per-source counts,
   and exit nonzero on any violation.

ID scheme:
  RMv2_NNNNNN   row number among the hom crossover rows of the cn then us
                Rodgers-Melnick tables, assigned before lift-over (gaps = dropped rows).
  EUROv2_NNNNNN row number in results/hmm_co_events_long.tsv (header excluded),
                assigned before lift-over.
  LRv4_NNNNNN / TEOv4_NNNNNN  row number in the lifted, sorted Samayoa v5 file
                (contiguous; sorted by chromosome number then start).

Needs results/hmm_co_events_long.tsv (scripts/hmm_co_pipeline.py) and CrossMap.
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).parent.parent
DATA = ROOT / "data"

RM_FILES = [DATA / "RodgersMelnick2015PNAS_cnnamImputedXOsegments.txt",
            DATA / "RodgersMelnick2015PNAS_usnamImputedXOsegments.txt"]
EURO_TSV = ROOT / "results/hmm_co_events_long.tsv"
SAMAYOA = {
    "LRv4": (DATA / "xo_Combined_LR13_LR14_parents2_AGPv4_filter0219.txt",
             "xo_Combined_LR13_LR14_parents2_v5.bed"),
    "TEOv4": (DATA / "xo_ZeaGBSv27raw_RareAllelesC2TeoCurated_depth_AGPv4_filtered0210.txt",
              "xo_ZeaGBSv27raw_RareAllelesC2TeoCurated_depth_v5.bed"),
}
V2V5_CHAIN = DATA / "v2v5.chain"
V4V5_CHAIN = DATA / "v4v5.chain"
V5_FAI = DATA / "Zm-B73-REFERENCE-NAM-5.0.fa.gz.fai"
CHROMS = [f"Chr{i}" for i in range(1, 11)]
SOURCES = ["RMv2", "EUROv2", "LRv4", "TEOv4"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Lift crossover sources to B73 v5 and build jri_v5.bed with validation.")
    parser.add_argument("--outdir", type=Path, default=DATA,
                        help="Directory for the four output BED files (default: data/).")
    parser.add_argument("--crossmap", default=shutil.which("CrossMap") or "/opt/anaconda3/bin/CrossMap",
                        help="CrossMap executable.")
    parser.add_argument("--workdir", type=Path,
                        help="Keep endpoint BED files here instead of a temporary directory.")
    return parser.parse_args()


def fail(msg):
    sys.exit(f"ERROR: {msg}")


def chrom_number(name):
    """'Chr7', 'chr7' or '7' -> 7; anything else (scaffolds) -> None."""
    m = re.fullmatch(r"(?:[Cc]hr)?(\d+)", name)
    n = int(m.group(1)) if m else None
    return n if n is not None and 1 <= n <= 10 else None


# ---------------------------------------------------------------- inputs

def read_rm():
    """Rodgers-Melnick hom crossovers, cn then us, as (chr, start, end, sample, id)."""
    rows = []
    for path in RM_FILES:
        for line in path.read_text().splitlines():
            line = line.replace("\r", "")
            if "het" in line or "Family" in line:  # drop header and het segments
                continue
            f = line.split("\t")
            rows.append((f"Chr{f[0]}", int(f[4]), int(f[5]), f[2], f"RMv2_{len(rows) + 1:06d}"))
    return rows


def read_euro():
    if not EURO_TSV.exists():
        fail(f"{EURO_TSV} missing; run scripts/hmm_co_pipeline.py --input-xlsx data/gb-2013-14-9-r103-S4.xlsx")
    lines = EURO_TSV.read_text().splitlines()
    header = lines[0].split("\t")
    col = {name: header.index(name) for name in
           ("sample_id", "chromosome", "left_coordinate", "right_coordinate")}
    rows = []
    for i, line in enumerate(lines[1:], start=1):
        f = line.split("\t")
        rows.append((f"Chr{f[col['chromosome']]}", int(f[col["left_coordinate"]]),
                     int(f[col["right_coordinate"]]), f[col["sample_id"]], f"EUROv2_{i:06d}"))
    return rows


def read_samayoa(path):
    """Samayoa rows as (chr, start, end, taxon, key); key = chr_start_end_taxon."""
    lines = path.read_text().splitlines()
    if lines[0].split("\t") != ["taxon", "chr", "start", "end"]:
        fail(f"unexpected header in {path}: {lines[0]!r}")
    rows = []
    for line in lines[1:]:
        taxon, chrom, start, end = line.split("\t")
        rows.append((chrom, int(start), int(end), taxon, f"{chrom}_{start}_{end}_{taxon}"))
    keys = Counter(r[4] for r in rows)
    dup = [k for k, n in keys.items() if n > 1]
    if dup:
        fail(f"{path.name}: {len(dup)} duplicated rows (would collapse during lift-over)")
    return rows


# ---------------------------------------------------------------- lift-over

def lift_endpoints(rows, chain, crossmap, workdir, tag):
    """Lift both endpoints of each interval; return {key: (v5 chr number, start, end)}.

    Endpoint markers are BED [start-1, start) and [end-1, end) for 1-based source
    positions. An interval is kept only if both markers lift, to the same
    chromosome 1..10; the result is a BED interval.
    """
    src = workdir / f"{tag}_markers_src.bed"
    dst = workdir / f"{tag}_markers_v5.bed"
    with open(src, "w") as fh:
        for chrom, start, end, _sample, key in rows:
            fh.write(f"{chrom}\t{start - 1}\t{start}\t{key}\n")
            fh.write(f"{chrom}\t{end - 1}\t{end}\t{key}\n")
    subprocess.run([crossmap, "bed", str(chain), str(src), str(dst)],
                   check=True, capture_output=True)

    total = Counter()
    hits = defaultdict(list)
    for line in dst.read_text().splitlines():
        chrom, start, end, key = line.split("\t")[:4]
        total[key] += 1
        hits[(key, chrom)].append((int(start), int(end)))

    lifted = {}
    for (key, chrom), ends in hits.items():
        n = chrom_number(chrom)
        if total[key] == 2 and len(ends) == 2 and n is not None:
            lifted[key] = (n, min(s for s, _ in ends), max(e for _, e in ends))
    return lifted


# ---------------------------------------------------------------- outputs

def fmt(row):
    return "\t".join(str(x) for x in row)


def write_lines(path, lines):
    path.write_text("".join(line + "\n" for line in lines))


def validate(path, lines, ncols, chrom_len, id_prefixes=None):
    """Check chromosome names, schema, coordinates and IDs of one output file."""
    seen_chroms = set()
    ids = Counter()
    for i, line in enumerate(lines, start=1):
        f = line.split("\t")
        if len(f) != ncols:
            fail(f"{path.name}:{i}: {len(f)} columns, expected {ncols}")
        chrom, start, end = f[0], int(f[1]), int(f[2])
        if chrom not in chrom_len:
            fail(f"{path.name}:{i}: chromosome {chrom!r} not in Chr1..Chr10")
        if not 0 <= start < end <= chrom_len[chrom]:
            fail(f"{path.name}:{i}: bad interval {chrom}:{start}-{end} (length {chrom_len[chrom]})")
        if not f[3]:
            fail(f"{path.name}:{i}: empty sample")
        seen_chroms.add(chrom)
        if id_prefixes:
            m = re.fullmatch(r"([A-Za-z0-9]+)_(\d{6})", f[4])
            if not m or m.group(1) not in id_prefixes:
                fail(f"{path.name}:{i}: bad ID {f[4]!r}")
            ids[f[4]] += 1
    if seen_chroms != set(CHROMS):
        fail(f"{path.name}: chromosomes {sorted(seen_chroms)} != Chr1..Chr10")
    dup = [k for k, n in ids.items() if n > 1]
    if dup:
        fail(f"{path.name}: {len(dup)} duplicate IDs, e.g. {dup[:3]}")
    return Counter(k.split("_")[0] for k in ids)


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    chrom_len = {}
    for line in V5_FAI.read_text().splitlines():
        name, length = line.split("\t")[:2]
        if name in {f"chr{i}" for i in range(1, 11)}:
            chrom_len["C" + name[1:]] = int(length)
    if set(chrom_len) != set(CHROMS):
        fail(f"{V5_FAI} lacks chr1..chr10")

    tmp = None
    if args.workdir:
        workdir = args.workdir
        workdir.mkdir(parents=True, exist_ok=True)
    else:
        tmp = tempfile.TemporaryDirectory()
        workdir = Path(tmp.name)

    raw = {"RMv2": read_rm(), "EUROv2": read_euro()}
    raw.update({src: read_samayoa(path) for src, (path, _) in SAMAYOA.items()})
    kept = {}

    # AGPv2 sources: IDs assigned before lift; sorted like `sort -k1,1 -k2,2n`.
    v2_rows = raw["RMv2"] + raw["EUROv2"]
    lifted = lift_endpoints(v2_rows, V2V5_CHAIN, args.crossmap, workdir, "v2")
    rm_euro = []
    for _chrom, _start, _end, sample, key in v2_rows:
        if key in lifted:
            n, s, e = lifted[key]
            rm_euro.append(fmt((f"Chr{n}", s, e, sample, key)))
    rm_euro.sort(key=lambda l: (l.split("\t", 1)[0].encode(), int(l.split("\t")[1]), l.encode()))
    out = args.outdir / "jri_rm_euro_v5.bed"
    write_lines(out, rm_euro)
    counts = validate(out, rm_euro, 5, chrom_len, {"RMv2", "EUROv2"})
    kept.update({src: counts[src] for src in ("RMv2", "EUROv2")})

    # AGPv4 sources: sort like `sort -k1,1n -k2,2n` on the bare-number
    # 5-column lines (chr, start, end, taxon, chr_start_end_taxon), then drop
    # the key, prefix Chr, and number rows in that order.
    combined = rm_euro[:]
    for src, (_path, name) in SAMAYOA.items():
        lifted = lift_endpoints(raw[src], V4V5_CHAIN, args.crossmap, workdir, src)
        rows = []
        for _chrom, _start, _end, taxon, key in raw[src]:
            if key in lifted:
                n, s, e = lifted[key]
                rows.append((n, s, fmt((n, s, e, taxon, key)).encode()))
        rows.sort()
        lines = []
        for n, _s, line in rows:
            f = line.decode().split("\t")
            lines.append(fmt((f"Chr{n}", f[1], f[2], f[3])))
        out = args.outdir / name
        write_lines(out, lines)
        validate(out, lines, 4, chrom_len)
        kept[src] = len(lines)
        combined += [f"{line}\t{src}_{i:06d}" for i, line in enumerate(lines, start=1)]

    # Step 3: `sort -k1,1V -k2,2n` (ties broken by whole-line byte order).
    combined.sort(key=lambda l: (int(l[3:l.index("\t")]), int(l.split("\t")[1]), l.encode()))
    out = args.outdir / "jri_v5.bed"
    write_lines(out, combined)
    counts = validate(out, combined, 5, chrom_len, set(SOURCES))

    print(f"{'source':<8} {'raw':>8} {'lifted':>8} {'in jri':>8} retention")
    for src in SOURCES:
        n_raw = len(raw[src])
        print(f"{src:<8} {n_raw:>8} {kept[src]:>8} {counts[src]:>8} {kept[src] / n_raw:8.1%}")
        if counts[src] != kept[src] or kept[src] == 0:
            fail(f"{src}: {counts[src]} IDs in jri_v5.bed vs {kept[src]} lifted intervals")
    print(f"total    {sum(len(r) for r in raw.values()):>8} {len(combined):>8}")
    print(f"wrote {args.outdir}/jri_rm_euro_v5.bed, xo_*_v5.bed and jri_v5.bed")
    if tmp:
        tmp.cleanup()


if __name__ == "__main__":
    main()
