#!/usr/bin/env python3
"""End-to-end release check for the tracked recombination-map products.

Verifies that the tracked products are mutually consistent and, with a
provenance manifest, that they were generated in dependency order from the
current inputs:

  HMM events -> lift-over (jri_v5.bed) -> finemap_v5 / hierarchical map
             -> HapMap exports / hotspot tracks

Consistency checks (always run):
  hmm           event count; every event's (sample_id, chromosome) is in the roster
  liftover      per-source audit counts and status categories; kept counts equal the
                ID-prefix counts in jri_v5.bed (and the per-source v5 files); the
                EUROv2 raw count equals the HMM event count
  jri           5 columns, Chr1..Chr10, 0 <= start < end <= length, unique IDs
  maps          6 columns, sorted non-overlapping segments, cM continuous and
                non-decreasing from 0, final cM equals the Ogut total (1e-6),
                cM_per_Mb consistent with the cM columns; finemap_v5 boundaries are
                jri_v5 interval endpoints
  hapmap        header, terminal rate 0, cM agrees with finemap_v5 at every segment
                boundary, last position equals the map end (or the chromosome end,
                when the builder pads to it at constant cM)
  hotspots      inside the finemap_v5 bounds

Provenance (content hashes, not timestamps):
  --write-manifest records sha256, size and row counts of every input, script and
  product plus the dependency graph. The default mode compares current files with
  the manifest and reports products that are stale, i.e. an upstream file changed
  since the manifest was written but the product did not.

Usage:
  python scripts/check_release.py                  # verify
  python scripts/check_release.py --write-manifest # record the current release
  python scripts/check_release.py --root DIR --manifest PATH
Exit status is nonzero when any check fails.
"""

import argparse
import csv
import datetime
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_ROOT = Path(__file__).resolve().parent.parent
CHROMS = [f"Chr{i}" for i in range(1, 11)]
SOURCES = ["RMv2", "EUROv2", "LRv4", "TEOv4"]
CM_TOL = 1e-6
MAX_LINES = 20

HMM_EVENTS = "results/hmm_co_events_long.tsv"
HMM_ROSTER = "results/hmm_sample_roster.tsv"
AUDIT = "results/liftover_audit.tsv"
JRI = "data/jri_v5.bed"
JRI_RM_EURO = "data/jri_rm_euro_v5.bed"
XO_V5 = {"LRv4": "data/xo_Combined_LR13_LR14_parents2_v5.bed",
         "TEOv4": "data/xo_ZeaGBSv27raw_RareAllelesC2TeoCurated_depth_v5.bed"}
FINEMAP = "data/finemap_v5.bed"
HIER = "data/finemap_hierarchical_v5.bed"
HAPMAP = [f"data/hapmap/chr{i}.hapmap.tsv" for i in range(1, 11)]
HOT30 = "data/hotspots_30x_v5.bed"
HOT1KB = "data/hotspots_1kb_20x_v5.bed"
OGUT = "data/ogut_fifthcM_map_agpv2.csv"
FAI = "data/v5.fa.gz.fai"
MANIFEST = "data/provenance_manifest.json"
HAPMAP_HEADER = "Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)"

# Dependency graph, in build order. Scripts count as inputs so code changes
# also mark their products stale.
GRAPH = {
    "hmm": {
        "inputs": ["data/gb-2013-14-9-r103-S4.xlsx", "scripts/hmm_co_pipeline.py"],
        "outputs": [HMM_EVENTS, HMM_ROSTER],
    },
    "liftover": {
        "inputs": [HMM_EVENTS,
                   "data/RodgersMelnick2015PNAS_cnnamImputedXOsegments.txt",
                   "data/RodgersMelnick2015PNAS_usnamImputedXOsegments.txt",
                   "data/xo_Combined_LR13_LR14_parents2_AGPv4_filter0219.txt",
                   "data/xo_ZeaGBSv27raw_RareAllelesC2TeoCurated_depth_AGPv4_filtered0210.txt",
                   "data/v2v5.chain", "data/v4v5.chain", FAI,
                   "scripts/build_jri_v5.py", "scripts/chain_liftover.py"],
        "outputs": [JRI, JRI_RM_EURO, XO_V5["LRv4"], XO_V5["TEOv4"], AUDIT],
    },
    "finemap": {
        "inputs": [JRI, OGUT, FAI, "scripts/build_finemap.py"],
        "outputs": [FINEMAP],
    },
    "hierarchical": {
        "inputs": [JRI, OGUT, FAI, "scripts/build_hierarchical_finemap.py"],
        "outputs": [HIER],
    },
    "hapmap": {
        "inputs": [FINEMAP, FAI, "scripts/build_finemap.py"],
        "outputs": HAPMAP,
    },
    "hotspots_30x": {
        "inputs": [FINEMAP, "scripts/define_hotspots.py"],
        "outputs": [HOT30],
    },
    "hotspots_1kb_20x": {
        "inputs": [FINEMAP, "scripts/define_hotspots_sliding.py"],
        "outputs": [HOT1KB],
    },
}


class CheckFailed(Exception):
    """Raised inside a check to abort it with a single failure message."""


class Report:
    def __init__(self, out=sys.stdout):
        self.out = out
        self.failed = []

    def add(self, name, failures, info="", warn=False):
        status = "PASS" if not failures else ("WARN" if warn else "FAIL")
        line = f"{status}  {name}"
        if info:
            line += f": {info}"
        print(line, file=self.out)
        for msg in failures[:MAX_LINES]:
            print(f"        - {msg}", file=self.out)
        if len(failures) > MAX_LINES:
            print(f"        - ... {len(failures) - MAX_LINES} more", file=self.out)
        if failures and not warn:
            self.failed.append(name)

    def run(self, name, func, *args):
        """Run func(*args) -> (failures, info); a CheckFailed or missing file fails it."""
        try:
            failures, info = func(*args)
        except CheckFailed as exc:
            failures, info = [str(exc)], ""
        except FileNotFoundError as exc:
            failures, info = [f"missing file: {exc.filename or exc}"], ""
        self.add(name, failures, info)
        return not failures


# ------------------------------------------------------------------ inputs

def require(root, rel):
    path = root / rel
    if not path.exists():
        raise CheckFailed(f"{rel} missing")
    return path


def chrom_lengths(root):
    lengths = {}
    for line in require(root, FAI).read_text().splitlines():
        name, length = line.split("\t")[:2]
        m = re.fullmatch(r"[Cc]hr(\d+)", name)
        if m and f"Chr{m.group(1)}" in CHROMS:
            lengths[f"Chr{m.group(1)}"] = int(length)
    if set(lengths) != set(CHROMS):
        raise CheckFailed(f"{FAI} lacks some of chr1..chr10")
    return lengths


def ogut_totals(root):
    """Per-chromosome Ogut cM span, derived exactly as build_finemap.py does."""
    ogut = pd.read_csv(require(root, OGUT))
    ogut["chr"] = "Chr" + ogut["chromosome"].astype(str)
    return {chrom: grp["cM"].max() - grp["cM"].min() for chrom, grp in ogut.groupby("chr")}


def read_jri(root):
    path = require(root, JRI)
    rows = [line.split("\t") for line in path.read_text().splitlines()]
    return rows


def source_of(audit_source, audit_id):
    """Audit source column is the lift tag ('v2' for RM+EURO); fall back to the ID prefix."""
    if audit_source in SOURCES:
        return audit_source
    prefix = audit_id.split("_", 1)[0]
    return prefix if prefix in SOURCES else audit_source


def count_lines(path):
    with open(path, "rb") as fh:
        return sum(1 for line in fh if line.strip())


# ------------------------------------------------------------------ checks

def check_hmm(root, ctx):
    events_path = require(root, HMM_EVENTS)
    roster_path = require(root, HMM_ROSTER)
    with open(roster_path, newline="") as fh:
        roster = list(csv.DictReader(fh, delimiter="\t"))
    if not roster or not {"sample_id", "chromosome"} <= set(roster[0]):
        raise CheckFailed(f"{HMM_ROSTER} lacks sample_id/chromosome columns")
    keys = Counter((r["sample_id"], r["chromosome"]) for r in roster)
    failures = []
    dup = [k for k, n in keys.items() if n > 1]
    if dup:
        failures.append(f"{len(dup)} duplicated roster rows, e.g. {dup[:3]}")
    with open(events_path, newline="") as fh:
        events = list(csv.DictReader(fh, delimiter="\t"))
    missing = [(e["sample_id"], e["chromosome"]) for e in events
               if (e["sample_id"], e["chromosome"]) not in keys]
    if missing:
        failures.append(f"{len(missing)} events whose sample/chromosome is not in the roster, "
                        f"e.g. {missing[:3]}")
    if not events:
        failures.append("no HMM events")
    ctx["hmm_events"] = len(events)
    ctx["hmm_rows"] = [(e["sample_id"], f"Chr{e['chromosome']}") for e in events]
    samples = len({r["sample_id"] for r in roster})
    return failures, (f"{len(events):,} events; roster {len(roster):,} sample-chromosomes, "
                      f"{samples:,} samples")


def check_liftover(root, ctx):
    path = require(root, AUDIT)
    raw, status = Counter(), defaultdict(Counter)
    ids = Counter()
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        col = header.index("status") if "status" in header else 8
        for row in reader:
            if not row:
                continue
            src = source_of(row[0], row[1])
            raw[src] += 1
            status[src][row[col]] += 1
            ids[(src, row[1])] += 1
    kept = {s: sum(n for k, n in status[s].items() if k.startswith("kept")) for s in SOURCES}
    failures = []
    dup = [k for k, n in ids.items() if n > 1]
    if dup:
        failures.append(f"{len(dup)} duplicated audit IDs, e.g. {dup[:3]}")
    unknown = sorted(set(raw) - set(SOURCES))
    if unknown:
        failures.append(f"unrecognised audit sources {unknown}")

    jri_counts = ctx.get("jri_prefix_counts")
    if jri_counts is None:
        jri_counts = Counter(r[4].split("_", 1)[0] for r in read_jri(root) if len(r) >= 5)
    for src in SOURCES:
        if raw[src] == 0:
            failures.append(f"{src}: no audit rows")
        if kept[src] != jri_counts.get(src, 0):
            failures.append(f"{src}: {kept[src]:,} kept in audit vs "
                            f"{jri_counts.get(src, 0):,} {src}_ IDs in {JRI}")
    for src, rel in XO_V5.items():
        if (root / rel).exists() and count_lines(root / rel) != kept[src]:
            failures.append(f"{rel}: {count_lines(root / rel):,} rows vs {kept[src]:,} kept {src}")
    if (root / JRI_RM_EURO).exists():
        n = count_lines(root / JRI_RM_EURO)
        if n != kept["RMv2"] + kept["EUROv2"]:
            failures.append(f"{JRI_RM_EURO}: {n:,} rows vs {kept['RMv2'] + kept['EUROv2']:,} "
                            "kept RMv2+EUROv2")
    if "hmm_events" in ctx and raw["EUROv2"] != ctx["hmm_events"]:
        failures.append(f"EUROv2: {raw['EUROv2']:,} audit rows vs {ctx['hmm_events']:,} "
                        f"events in {HMM_EVENTS}")

    lines = []
    for src in SOURCES:
        cats = ", ".join(f"{k} {n:,}" for k, n in sorted(status[src].items(), key=lambda x: -x[1]))
        lines.append(f"{src} raw {raw[src]:,} kept {kept[src]:,} ({cats})")
    return failures, "\n          " + "\n          ".join(lines)


def check_jri(root, ctx):
    lengths = chrom_lengths(root)
    rows = read_jri(root)
    failures = []
    ids = Counter()
    seen = set()
    for i, f in enumerate(rows, start=1):
        if len(f) != 5:
            failures.append(f"line {i}: {len(f)} columns, expected 5")
            continue
        chrom = f[0]
        if chrom not in lengths:
            failures.append(f"line {i}: chromosome {chrom!r} not in Chr1..Chr10")
            continue
        start, end = int(f[1]), int(f[2])
        if not 0 <= start < end <= lengths[chrom]:
            failures.append(f"line {i}: bad interval {chrom}:{start}-{end} "
                            f"(length {lengths[chrom]})")
        if not re.fullmatch(r"(RMv2|EUROv2|LRv4|TEOv4)_\d{6}", f[4]):
            failures.append(f"line {i}: bad ID {f[4]!r}")
        ids[f[4]] += 1
        seen.add(chrom)
    dup = [k for k, n in ids.items() if n > 1]
    if dup:
        failures.append(f"{len(dup)} duplicate IDs, e.g. {dup[:3]}")
    if seen != set(CHROMS) and not failures:
        failures.append(f"chromosomes present {sorted(seen)} != Chr1..Chr10")
    prefix = Counter(k.split("_", 1)[0] for k in ids)
    ctx["jri_prefix_counts"] = prefix

    # EUROv2_NNNNNN is row NNNNNN of the HMM events table (header excluded); a
    # jri_v5.bed built from another HMM run has mismatched samples/chromosomes.
    events = ctx.get("hmm_rows")
    if events is not None:
        bad = []
        for f in rows:
            if len(f) == 5 and f[4].startswith("EUROv2_"):
                n = int(f[4].split("_")[1])
                if n > len(events) or events[n - 1] != (f[3], f[0]):
                    bad.append(f[4])
        if bad:
            failures.append(f"{len(bad):,} of {prefix.get('EUROv2', 0):,} EUROv2 intervals do not "
                            f"match their row in {HMM_EVENTS} (sample/chromosome), e.g. {bad[:3]}; "
                            "jri_v5.bed was built from a different HMM run")
    ends = defaultdict(set)
    for f in rows:
        if len(f) == 5:
            ends[f[0]].update((int(f[1]), int(f[2])))
    ctx["jri_endpoints"] = ends
    return failures, (f"{len(rows):,} intervals (" +
                      ", ".join(f"{s} {prefix.get(s, 0):,}" for s in SOURCES) + ")")


def load_map(root, rel):
    path = require(root, rel)
    df = pd.read_csv(path, sep="\t", header=None)
    if df.shape[1] != 6:
        raise CheckFailed(f"{rel}: {df.shape[1]} columns, expected 6")
    df.columns = ["chr", "start", "end", "cm_start", "cm_end", "rate"]
    return df


def check_map(root, ctx, rel):
    lengths = chrom_lengths(root)
    totals = ogut_totals(root)
    df = load_map(root, rel)
    failures = []
    bad_chr = sorted(set(df["chr"]) - set(CHROMS))
    if bad_chr:
        raise CheckFailed(f"chromosomes outside Chr1..Chr10: {bad_chr}")
    if set(df["chr"]) != set(CHROMS):
        failures.append(f"missing chromosomes {sorted(set(CHROMS) - set(df['chr']))}")
    blocks = (df["chr"] != df["chr"].shift()).sum()
    if blocks != df["chr"].nunique():
        failures.append("chromosome rows are not contiguous")

    bounds, final = {}, {}
    for chrom, g in df.groupby("chr", sort=False):
        s, e = g["start"].to_numpy(), g["end"].to_numpy()
        cs, ce, r = (g[c].to_numpy(dtype=float) for c in ("cm_start", "cm_end", "rate"))
        if (s < 0).any() or (e > lengths[chrom]).any():
            failures.append(f"{chrom}: segment outside 0..{lengths[chrom]}")
        if (e <= s).any():
            failures.append(f"{chrom}: {(e <= s).sum()} empty or reversed segments")
        if (s[1:] < e[:-1]).any():
            i = int(np.argmax(s[1:] < e[:-1]))
            failures.append(f"{chrom}: unsorted/overlapping segments at {s[i + 1]} < {e[i]}")
        if abs(cs[0]) > 1e-12:
            failures.append(f"{chrom}: first cM_start {cs[0]} != 0")
        if (ce < cs).any():
            i = int(np.argmax(ce < cs))
            failures.append(f"{chrom}: decreasing cM at {s[i]}-{e[i]} ({cs[i]} -> {ce[i]})")
        gap = np.abs(cs[1:] - ce[:-1])
        if (gap > 1e-9 * np.maximum(1.0, np.abs(ce[:-1]))).any():
            i = int(np.argmax(gap))
            failures.append(f"{chrom}: cM discontinuity at {s[i + 1]} "
                            f"({ce[i]} -> {cs[i + 1]})")
        implied = r * (e - s) / 1e6
        dcm = ce - cs
        bad = np.abs(implied - dcm) > 1e-9 + 1e-6 * np.abs(dcm)
        if bad.any():
            i = int(np.argmax(bad))
            failures.append(f"{chrom}: {bad.sum()} segments with cM_per_Mb inconsistent with cM, "
                            f"e.g. {s[i]}-{e[i]} rate {r[i]} vs {dcm[i] / (e[i] - s[i]) * 1e6}")
        target = totals.get(chrom)
        if target is None:
            failures.append(f"{chrom}: no Ogut total")
        elif abs(ce[-1] - target) > CM_TOL:
            failures.append(f"{chrom}: final cM {ce[-1]:.9f} != Ogut total {target:.9f}")
        bounds[chrom] = (int(s[0]), int(e[-1]))
        final[chrom] = float(ce[-1])

    if rel == FINEMAP:
        ctx["finemap"] = df
        ctx["finemap_bounds"] = bounds
        ends = ctx.get("jri_endpoints")
        if ends:
            for chrom, g in df.groupby("chr", sort=False):
                pos = set(g["start"]).union(g["end"])
                foreign = pos - ends.get(chrom, set())
                if foreign:
                    failures.append(f"{chrom}: {len(foreign)} segment boundaries are not "
                                    f"{JRI} endpoints, e.g. {sorted(foreign)[:3]}")
                jri_lo, jri_hi = min(ends[chrom]), max(ends[chrom])
                if bounds.get(chrom) != (jri_lo, jri_hi):
                    failures.append(f"{chrom}: map spans {bounds.get(chrom)} but {JRI} spans "
                                    f"{(jri_lo, jri_hi)}")
    return failures, f"{len(df):,} segments"


def check_hapmap(root, ctx):
    lengths = chrom_lengths(root)
    fm = ctx.get("finemap")
    if fm is None:
        fm = load_map(root, FINEMAP)
    failures = []
    total_rows = 0
    for i, rel in enumerate(HAPMAP, start=1):
        chrom = f"Chr{i}"
        path = root / rel
        if not path.exists():
            failures.append(f"{rel} missing")
            continue
        with open(path) as fh:
            header = fh.readline().rstrip("\n")
        if header != HAPMAP_HEADER:
            failures.append(f"{rel}: bad header {header!r}")
            continue
        hm = pd.read_csv(path, sep="\t")
        hm.columns = ["chrom", "pos", "rate", "cm"]
        total_rows += len(hm)
        pos, cm = hm["pos"].to_numpy(), hm["cm"].to_numpy(dtype=float)
        if (hm["chrom"] != f"chr{i}").any():
            failures.append(f"{rel}: chromosome column is not chr{i}")
        if (np.diff(pos) <= 0).any():
            failures.append(f"{rel}: positions not strictly increasing")
        if (np.diff(cm) < -1e-9).any():
            failures.append(f"{rel}: cM decreases")
        if hm["rate"].iloc[-1] != 0:
            failures.append(f"{rel}: terminal rate {hm['rate'].iloc[-1]} != 0")
        seg = fm[fm["chr"] == chrom]
        if seg.empty:
            failures.append(f"{rel}: no {chrom} segments in {FINEMAP}")
            continue
        bpos = np.concatenate([seg["start"].to_numpy(), seg["end"].to_numpy()])
        bcm = np.concatenate([seg["cm_start"].to_numpy(float), seg["cm_end"].to_numpy(float)])
        idx = np.clip(np.searchsorted(pos, bpos), 0, len(pos) - 1)
        absent = pos[idx] != bpos
        if absent.any():
            failures.append(f"{rel}: {absent.sum()} map boundaries absent, "
                            f"e.g. {bpos[absent][:3].tolist()}")
        diff = np.abs(cm[idx] - bcm)
        off = ~absent & (diff > CM_TOL + 1e-9 * np.abs(bcm))
        if off.any():
            j = int(np.argmax(np.where(off, diff, -1)))
            failures.append(f"{rel}: {off.sum()} boundaries disagree with {FINEMAP}, "
                            f"e.g. {bpos[j]}: {cm[idx][j]} vs {bcm[j]}")
        map_end, map_cm = int(seg["end"].iloc[-1]), float(seg["cm_end"].iloc[-1])
        if pos[-1] not in (map_end, lengths[chrom]):
            failures.append(f"{rel}: last position {pos[-1]} != map end {map_end} "
                            f"(or chromosome end {lengths[chrom]})")
        if abs(cm[-1] - map_cm) > CM_TOL + 1e-9 * abs(map_cm):
            failures.append(f"{rel}: last cM {cm[-1]} != map final cM {map_cm}")
    return failures, f"{len(HAPMAP)} files, {total_rows:,} rows"


def check_hotspots(root, ctx, rel):
    bounds = ctx.get("finemap_bounds")
    if bounds is None:
        df = load_map(root, FINEMAP)
        bounds = {c: (int(g["start"].iloc[0]), int(g["end"].iloc[-1]))
                  for c, g in df.groupby("chr", sort=False)}
    failures = []
    n = 0
    for i, line in enumerate(require(root, rel).read_text().splitlines(), start=1):
        f = line.split("\t")
        n += 1
        chrom = "C" + f[0][1:] if f[0].startswith("chr") else f[0]
        if chrom not in bounds:
            failures.append(f"line {i}: chromosome {f[0]!r} not in the map")
            continue
        s, e = int(f[1]), int(f[2])
        lo, hi = bounds[chrom]
        if not lo <= s < e <= hi:
            failures.append(f"line {i}: {f[0]}:{s}-{e} outside map bounds {lo}-{hi}")
    return failures, f"{n:,} regions"


# ------------------------------------------------------------------ manifest

def file_record(root, rel):
    path = root / rel
    if not path.exists():
        return None
    digest = hashlib.sha256()
    rows = 0
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
            rows += chunk.count(b"\n")
    text = path.suffix not in {".xlsx", ".gz", ".npz", ".png"}
    return {"sha256": digest.hexdigest(), "bytes": path.stat().st_size,
            "rows": rows if text else None}


def all_files(graph):
    seen = []
    for node in graph.values():
        for rel in node["inputs"] + node["outputs"]:
            if rel not in seen:
                seen.append(rel)
    return seen


def ancestors(graph, rel):
    """All files upstream of rel (transitively)."""
    producers = {out: node for node in graph.values() for out in node["outputs"]}
    result, stack = set(), [rel]
    while stack:
        node = producers.get(stack.pop())
        if not node:
            continue
        for inp in node["inputs"]:
            if inp not in result:
                result.add(inp)
                stack.append(inp)
    return result


def write_manifest(root, manifest_path):
    records = {rel: file_record(root, rel) for rel in all_files(GRAPH)}
    manifest = {
        "version": 1,
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "graph": GRAPH,
        "files": records,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=1) + "\n")
    return records


def check_manifest(root, manifest_path, report):
    if not manifest_path.exists():
        report.add("provenance manifest", ["no manifest; run with --write-manifest once the "
                                           "consistency checks pass"],
                   f"{manifest_path} not found, hash/dependency checks skipped", warn=True)
        return
    manifest = json.loads(manifest_path.read_text())
    graph = manifest.get("graph", GRAPH)
    recorded = manifest.get("files", {})
    changed = []
    for rel in all_files(graph):
        old = recorded.get(rel)
        new = file_record(root, rel)
        if (old or {}).get("sha256") != (new or {}).get("sha256"):
            changed.append(rel)
    products = [rel for node in graph.values() for rel in node["outputs"]]
    stale, rebuilt = [], []
    for rel in products:
        up = sorted(a for a in ancestors(graph, rel) if a in changed)
        if up and rel not in changed:
            stale.append(f"{rel} is stale (upstream changed: {', '.join(up)})")
        elif rel in changed and not up:
            rebuilt.append(f"{rel} changed although none of its upstream files did")
    changed_msgs = [f"{rel} differs from the manifest" for rel in changed]
    report.add("provenance manifest", changed_msgs,
               f"{manifest_path.name} written {manifest.get('written', '?')}, "
               f"{len(recorded)} files, {len(changed)} changed")
    report.add("dependency order", stale + rebuilt,
               "no stale products" if not (stale or rebuilt)
               else f"{len(stale)} stale products")


# ------------------------------------------------------------------ main

def run_consistency(root, report):
    ctx = {}
    report.run("hmm events/roster", check_hmm, root, ctx)
    report.run("jri_v5.bed", check_jri, root, ctx)
    report.run("lift-over audit", check_liftover, root, ctx)
    report.run("finemap_v5.bed", check_map, root, ctx, FINEMAP)
    report.run("finemap_hierarchical_v5.bed", check_map, root, ctx, HIER)
    report.run("hapmap exports", check_hapmap, root, ctx)
    report.run("hotspots_30x_v5.bed", check_hotspots, root, ctx, HOT30)
    report.run("hotspots_1kb_20x_v5.bed", check_hotspots, root, ctx, HOT1KB)


def main(argv=None, out=sys.stdout):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                        help="repository root (default: this repository)")
    parser.add_argument("--manifest", type=Path,
                        help=f"manifest path (default: ROOT/{MANIFEST})")
    parser.add_argument("--write-manifest", action="store_true",
                        help="record hashes of the current inputs and products")
    parser.add_argument("--force", action="store_true",
                        help="write the manifest even if consistency checks fail")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    manifest_path = args.manifest or root / MANIFEST

    report = Report(out)
    run_consistency(root, report)
    if args.write_manifest:
        if report.failed and not args.force:
            print(f"\nFAIL  not writing {manifest_path}: {len(report.failed)} consistency "
                  "checks failed (use --force to override)", file=out)
            return 1
        records = write_manifest(root, manifest_path)
        missing = [rel for rel, rec in records.items() if rec is None]
        print(f"\nwrote {manifest_path} ({len(records)} files"
              + (f"; absent: {', '.join(missing)}" if missing else "") + ")", file=out)
    else:
        check_manifest(root, manifest_path, report)

    if report.failed:
        print(f"\nRELEASE CHECK FAILED ({len(report.failed)}): {', '.join(report.failed)}",
              file=out)
        return 1
    print("\nRELEASE CHECK PASSED", file=out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
