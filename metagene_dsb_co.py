#!/usr/bin/env python3

#### module load python/3.9.19

#### USAGE:
####   script.py --gff genes.gff3 --bed1 co.bed --bed2 dsb.bed --flank 10000 20 --gene 2000 20
####   script.py --gene_bed genes.bed --bed1 co.bed --flank 10000 20 --gene 2000 20

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import sys
from pathlib import Path

# =====================================================
# usage system & input files
# =====================================================


def parse_args():
    parser = argparse.ArgumentParser(description="Plot dual-anchor metagene profiles from one or two BED signal tracks.")
    gene_group = parser.add_mutually_exclusive_group()
    gene_group.add_argument("--gff", dest="gff_file", help="Gene annotation in GFF/GFF3 format.")
    gene_group.add_argument("--gene_bed", dest="gene_bed_file", help="Gene annotation in the original BED-like format.")
    parser.add_argument("--bed1", required=True, help="First signal BED with columns: chr, start, end, count.")
    parser.add_argument("--bed2", help="Optional second signal BED with columns: chr, start, end, count.")
    parser.add_argument("--flank", type=int, nargs=2, metavar=("SIZE_BP", "N_BINS"), required=True,
                        help="Flank region size in bp and number of bins.")

    parser.add_argument("--gene", type=int, nargs=2, metavar=("SIZE_BP", "N_BINS"), required=True,
                        help="Inner gene-boundary region size in bp and number of bins.")
    parser.add_argument("--uniform", action="store_true",
                        help="Treat each BED interval as one event spread uniformly across its span instead of using midpoint counts.")
    return parser.parse_args()


args = parse_args()

GENE_FILE = args.gff_file or args.gene_bed_file
if GENE_FILE is None:
    raise SystemExit("You must provide either --gff or --gene_bed")
CO_FILE = args.bed1
DSB_FILE = args.bed2
FLANK_BIN = args.flank
GENE_BIN = args.gene
UNIFORM = args.uniform


_warned_chrom_renames = set()


def normalize_chrom_name(value):
    if pd.isna(value):
        return value
    original = str(value).strip()
    lower = original.lower()
    if lower.startswith("chr"):
        normalized = f"chr{original[3:]}"
    else:
        normalized = original
    if normalized != original and original not in _warned_chrom_renames:
        print(
            f"Warning: normalized chromosome name '{original}' -> '{normalized}'",
            file=sys.stderr,
        )
        _warned_chrom_renames.add(original)
    return normalized

# region sizes
FLANK_SIZE, FLANK_N_BINS = FLANK_BIN
INNER_SIZE, INNER_N_BINS = GENE_BIN

SMOOTH = 5

# =====================================================
# define metagene layout
# =====================================================

n_up    = FLANK_N_BINS
n_inner = INNER_N_BINS
n_exon  = 1

TOTAL = n_up + n_inner + n_exon + n_exon + n_inner + n_up

# offsets
OFF_TSS  = n_up
OFF_FIRST_EXON = OFF_TSS + n_inner
OFF_LAST_EXON = OFF_FIRST_EXON + n_exon
OFF_TTS  = OFF_LAST_EXON + n_exon
OFF_DOWN = OFF_TTS + n_inner

# =====================================================
# load & process data
# =====================================================

if GENE_FILE.endswith((".gff3", ".gff", ".gff3.gz", ".gff.gz")):
    genes = pd.read_csv(
        GENE_FILE,
        sep="\t",
        header=None,
        comment="#",
        names=["chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"],
        usecols=list(range(9)),
    )
    genes = genes[genes["feature"] == "gene"].copy()
    genes["start"] = genes["start"] - 1
    genes["len"] = genes["end"] - genes["start"]
    def parse_attrs(attr_str):
        result = {}
        for item in str(attr_str).split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                result[key] = value
        return result

    genes["attr_map"] = genes["attributes"].map(parse_attrs)
    genes["gene_id"] = genes["attr_map"].map(lambda x: x.get("ID"))
    genes["canonical_transcript"] = genes["attr_map"].map(lambda x: x.get("canonical_transcript"))

    exons = pd.read_csv(
        GENE_FILE,
        sep="\t",
        header=None,
        comment="#",
        names=["chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"],
        usecols=list(range(9)),
    )
    exons = exons[exons["feature"] == "exon"].copy()
    exons["start"] = exons["start"] - 1
    exons["attr_map"] = exons["attributes"].map(parse_attrs)
    exons["parent"] = exons["attr_map"].map(lambda x: x.get("Parent"))
    exons["rank"] = pd.to_numeric(exons["attr_map"].map(lambda x: x.get("rank")), errors="coerce")
    exons["gene_id"] = exons["parent"].map(
        lambda x: x.rsplit("_T", 1)[0] if isinstance(x, str) and "_T" in x else None
    )

    exon_summary = {}
    for gene_id, sub in exons.dropna(subset=["gene_id"]).groupby("gene_id", sort=False):
        transcript_groups = {tid: grp.copy() for tid, grp in sub.groupby("parent", sort=False)}
        canonical = genes.loc[genes["gene_id"] == gene_id, "canonical_transcript"]
        canonical_id = canonical.iloc[0] if not canonical.empty else None
        if canonical_id in transcript_groups:
            chosen = transcript_groups[canonical_id]
        else:
            chosen = max(transcript_groups.values(), key=lambda grp: len(grp))
        if chosen["rank"].notna().any():
            chosen = chosen.sort_values("rank", kind="mergesort")
        else:
            chosen = chosen.sort_values("start", kind="mergesort")
        first_row = chosen.iloc[0]
        last_row = chosen.iloc[-1]
        exon_summary[gene_id] = (
            int(first_row["start"]),
            int(first_row["end"]),
            int(last_row["start"]),
            int(last_row["end"]),
        )

    genes["first_exon_start"] = genes["gene_id"].map(lambda gid: exon_summary.get(gid, (np.nan,)*4)[0])
    genes["first_exon_end"] = genes["gene_id"].map(lambda gid: exon_summary.get(gid, (np.nan,)*4)[1])
    genes["last_exon_start"] = genes["gene_id"].map(lambda gid: exon_summary.get(gid, (np.nan,)*4)[2])
    genes["last_exon_end"] = genes["gene_id"].map(lambda gid: exon_summary.get(gid, (np.nan,)*4)[3])
    genes = genes.dropna(subset=["first_exon_start", "first_exon_end", "last_exon_start", "last_exon_end"]).copy()
    genes["first_exon_start"] = genes["first_exon_start"].astype(np.int64)
    genes["first_exon_end"] = genes["first_exon_end"].astype(np.int64)
    genes["last_exon_start"] = genes["last_exon_start"].astype(np.int64)
    genes["last_exon_end"] = genes["last_exon_end"].astype(np.int64)
else:
    genes = pd.read_csv(
        GENE_FILE, sep="\t", header=None,
        names=["chr","start","end","len","dot","strand","type"]
    )

def load_signal_df(path):
    if path is None:
        return None
    df = pd.read_csv(path, sep="	", header=None)
    if df.shape[1] < 3:
        raise SystemExit(f"Signal BED must have at least 3 columns: {path}")
    df = df.iloc[:, :4].copy()
    cols = ["chr", "start", "end"]
    if df.shape[1] >= 4:
        cols.append("count")
    df.columns = cols
    if "count" not in df.columns:
        df["count"] = 1.0
    return df


co = load_signal_df(CO_FILE)
dsb = load_signal_df(DSB_FILE)


def normalize_signal_df(df):
    if df is None:
        return None
    df = df.copy()
    df["chr"] = df["chr"].map(normalize_chrom_name)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(1.0)
    df = df.dropna(subset=["chr", "start", "end"])
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    return df


co = normalize_signal_df(co)
dsb = normalize_signal_df(dsb)

if not UNIFORM:
    co["mid"] = ((co.start + co.end) // 2).astype(np.int64)
    if dsb is not None:
        dsb["mid"] = ((dsb.start + dsb.end) // 2).astype(np.int64)


def build_signal_dict(df):
    signal = {}
    for chrom, sub in df.groupby("chr", sort=False):
        if UNIFORM:
            sub = sub.sort_values("start", kind="mergesort").copy()
            lengths = (sub["end"] - sub["start"]).to_numpy(dtype=np.int64)
            valid = lengths > 0
            signal[chrom] = (
                sub.loc[valid, "start"].to_numpy(dtype=np.int64),
                sub.loc[valid, "end"].to_numpy(dtype=np.int64),
                1.0 / lengths[valid],
            )
        else:
            sub = sub.sort_values("mid", kind="mergesort")
            signal[chrom] = (
                sub["mid"].to_numpy(dtype=np.int64),
                sub["count"].to_numpy(),
            )
    return signal


co_chr = build_signal_dict(co)
dsb_chr = build_signal_dict(dsb) if dsb is not None else None

# =====================================================
# gene profile function
# =====================================================

def gene_profile(g, signal_dict):

    left = np.zeros(OFF_TTS)
    right = np.zeros(TOTAL - OFF_TTS)
    left_widths = np.zeros(OFF_TTS)
    right_widths = np.zeros(TOTAL - OFF_TTS)
    if g.chr not in signal_dict:
        return left, 0.0, 0.0, right

    if UNIFORM:
        starts, ends, weights = signal_dict[g.chr]

        def range_signal(lo, hi):
            left_idx = max(0, np.searchsorted(starts, lo, side="right") - 1)
            right_idx = np.searchsorted(starts, hi, side="left")
            total = 0.0
            for idx in range(left_idx, right_idx):
                ov_start = max(lo, starts[idx])
                ov_end = min(hi, ends[idx])
                if ov_end > ov_start:
                    total += (ov_end - ov_start) * weights[idx]
            return total
    else:
        m, c = signal_dict[g.chr]

        def slice_range(lo, hi):
            left_idx = np.searchsorted(m, lo, side="left")
            right_idx = np.searchsorted(m, hi, side="left")
            return m[left_idx:right_idx], c[left_idx:right_idx]

    if g.strand == "+":
        TSS, TTS = g.start, g.end
        first_exon = (g.first_exon_start, g.first_exon_end)
        last_exon = (g.last_exon_start, g.last_exon_end)
    else:
        TSS, TTS = g.end, g.start
        first_exon = (g.last_exon_start, g.last_exon_end)
        last_exon = (g.first_exon_start, g.first_exon_end)

    # upstream flank
    if UNIFORM:
        for i in range(n_up):
            b0 = TSS - FLANK_SIZE + i * FLANK_SIZE / n_up
            b1 = TSS - FLANK_SIZE + (i + 1) * FLANK_SIZE / n_up
            lo = int(round(b0))
            hi = int(round(b1))
            left[i] += range_signal(lo, hi)
            left_widths[i] += max(0, hi - lo)
    else:
        m_sub, c_sub = slice_range(TSS - FLANK_SIZE, TSS)
        idx = (((m_sub - (TSS - FLANK_SIZE)) / FLANK_SIZE) * n_up).astype(int)
        valid = (idx >= 0) & (idx < n_up)
        np.add.at(left, idx[valid], c_sub[valid])
        left_widths[:n_up] += FLANK_SIZE / n_up

    # internal TSS window
    if UNIFORM:
        for i in range(n_inner):
            b0 = TSS + i * INNER_SIZE / n_inner
            b1 = TSS + (i + 1) * INNER_SIZE / n_inner
            lo = int(round(b0))
            hi = int(round(b1))
            left[OFF_TSS + i] += range_signal(lo, hi)
            left_widths[OFF_TSS + i] += max(0, hi - lo)
    else:
        m_sub, c_sub = slice_range(TSS, TSS + INNER_SIZE)
        idx = OFF_TSS + (((m_sub - TSS) / INNER_SIZE) * n_inner).astype(int)
        valid = (idx >= OFF_TSS) & (idx < OFF_TTS)
        np.add.at(left, idx[valid], c_sub[valid])
        left_widths[OFF_TSS:OFF_TTS] += INNER_SIZE / n_inner

    # first exon point
    exon_start, exon_end = first_exon
    exon_width = max(1, exon_end - exon_start)
    if UNIFORM:
        first_exon_value = range_signal(exon_start, exon_end) / exon_width
    else:
        _m_sub, c_sub = slice_range(exon_start, exon_end)
        first_exon_value = c_sub.sum() / exon_width

    # last exon point
    exon_start, exon_end = last_exon
    exon_width = max(1, exon_end - exon_start)
    if UNIFORM:
        last_exon_value = range_signal(exon_start, exon_end) / exon_width
    else:
        _m_sub, c_sub = slice_range(exon_start, exon_end)
        last_exon_value = c_sub.sum() / exon_width

    # internal TTS window
    if UNIFORM:
        for i in range(n_inner):
            b0 = TTS - INNER_SIZE + i * INNER_SIZE / n_inner
            b1 = TTS - INNER_SIZE + (i + 1) * INNER_SIZE / n_inner
            lo = int(round(b0))
            hi = int(round(b1))
            right[i] += range_signal(lo, hi)
            right_widths[i] += max(0, hi - lo)
    else:
        m_sub, c_sub = slice_range(TTS - INNER_SIZE, TTS)
        idx = (((m_sub - (TTS - INNER_SIZE)) / INNER_SIZE) * n_inner).astype(int)
        valid = (idx >= 0) & (idx < n_inner)
        np.add.at(right, idx[valid], c_sub[valid])
        right_widths[:n_inner] += INNER_SIZE / n_inner

    # downstream flank
    if UNIFORM:
        for i in range(n_up):
            b0 = TTS + i * FLANK_SIZE / n_up
            b1 = TTS + (i + 1) * FLANK_SIZE / n_up
            lo = int(round(b0))
            hi = int(round(b1))
            right[n_inner + i] += range_signal(lo, hi)
            right_widths[n_inner + i] += max(0, hi - lo)
    else:
        m_sub, c_sub = slice_range(TTS, TTS + FLANK_SIZE)
        idx = n_inner + (((m_sub - TTS) / FLANK_SIZE) * n_up).astype(int)
        valid = (idx >= n_inner) & (idx < len(right))
        np.add.at(right, idx[valid], c_sub[valid])
        right_widths[n_inner:] += FLANK_SIZE / n_up

    left = np.divide(left, left_widths, out=np.zeros_like(left), where=left_widths > 0)
    right = np.divide(right, right_widths, out=np.zeros_like(right), where=right_widths > 0)
    return left, first_exon_value, last_exon_value, right

# =====================================================
# aggregation function
# =====================================================

def aggregate(signal_dict):

    left_profiles = []
    first_exon_values = []
    last_exon_values = []
    right_profiles = []
    excluded_short = 0

    for g in genes.itertuples(index=False):
        if g.len <= 2 * INNER_SIZE:
            excluded_short += 1
            continue
        left, first_exon_value, last_exon_value, right = gene_profile(g, signal_dict)
        if left.sum() == 0 and right.sum() == 0 and first_exon_value == 0 and last_exon_value == 0:
            continue
        left_profiles.append(left)
        first_exon_values.append(first_exon_value)
        last_exon_values.append(last_exon_value)
        right_profiles.append(right)

    print(
        f"{excluded_short} genes excluded because they are shorter than 2*INNER_SIZE",
        file=sys.stderr,
    )

    if not left_profiles:
        raise SystemExit(
            "No overlapping signal was found between the gene annotations and signal BED files after chromosome normalization."
        )

    return (
        np.vstack(left_profiles).mean(axis=0),
        float(np.mean(first_exon_values)),
        float(np.mean(last_exon_values)),
        np.vstack(right_profiles).mean(axis=0),
    )

# =====================================================
# smoothing function
# =====================================================

def smooth(x, w=3):
    return np.convolve(x, np.ones(w)/w, mode="same")

# =====================================================
# run aggregation & smoothing
# =====================================================

co_left_raw, co_first_exon, co_last_exon, co_right_raw = aggregate(co_chr)
dsb_parts = aggregate(dsb_chr) if dsb_chr is not None else None
if dsb_parts is not None:
    dsb_left_raw, dsb_first_exon, dsb_last_exon, dsb_right_raw = dsb_parts

co_left_raw = smooth(co_left_raw, SMOOTH)
co_right_raw = smooth(co_right_raw, SMOOTH)
if dsb_parts is not None:
    dsb_left_raw = smooth(dsb_left_raw, SMOOTH)
    dsb_right_raw = smooth(dsb_right_raw, SMOOTH)

# -----------------------------------------------------
# plot: define boundaries
# -----------------------------------------------------

LEFT_TOTAL = OFF_TTS
RIGHT_TOTAL = n_inner + n_up

# -----------------------------------------------------
# plot: helper function to format x-axis
# -----------------------------------------------------

def ticks_for_region(start_pos, end_pos, n_bins, max_ticks=5):
    if n_bins <= 1:
        return [0, n_bins - 1], [start_pos, end_pos]
    step = max(1, int(np.ceil((n_bins - 1) / (max_ticks - 1))))
    idxs = list(range(0, n_bins, step))
    if idxs[-1] != n_bins - 1:
        idxs.append(n_bins - 1)
    labels = []
    for idx in idxs:
        frac = idx / (n_bins - 1) if n_bins > 1 else 0
        labels.append(start_pos + frac * (end_pos - start_pos))
    return idxs, labels


def fmt_bp(value_bp):
    if abs(value_bp) >= 1000:
        kb = value_bp / 1000
        if abs(kb - round(kb)) < 1e-9:
            return f"{int(round(kb))} kb"
        return f"{kb:g} kb"
    return f"{int(round(value_bp))} bp"


def style_dual_anchor(ax_left, ax_mid, ax_right):
    up_ticks, up_labels = ticks_for_region(-FLANK_SIZE, 0, n_up, max_ticks=4)
    tss_ticks, tss_labels = ticks_for_region(0, INNER_SIZE, n_inner, max_ticks=3)
    left_ticks = up_ticks[:-1] + [OFF_TSS + t for t in tss_ticks]
    left_labels = [fmt_bp(v) for v in up_labels[:-1]] + ["TSS"] + [f"+{fmt_bp(v)}" for v in tss_labels[1:]]
    ax_left.set_xticks(left_ticks)
    ax_left.set_xticklabels(left_labels, rotation=45, ha="right")
    ax_left.axvspan(OFF_TSS, LEFT_TOTAL, alpha=0.08)
    ax_left.axvline(OFF_TSS, linestyle="--", linewidth=1)
    ax_left.set_xlim(0, LEFT_TOTAL - 1)

    ax_mid.set_xticks([0, 1])
    ax_mid.set_xticklabels(["First exon", "Last exon"], rotation=45, ha="right")
    ax_mid.set_xlim(-0.5, 1.5)

    tts_ticks, tts_labels = ticks_for_region(-INNER_SIZE, 0, n_inner, max_ticks=3)
    down_ticks, down_labels = ticks_for_region(0, FLANK_SIZE, n_up, max_ticks=4)
    right_ticks = tts_ticks[:-1] + [n_inner + t for t in down_ticks]
    right_labels = [f"-{fmt_bp(abs(v))}" for v in tts_labels[:-1]] + ["TTS"] + [f"+{fmt_bp(v)}" for v in down_labels[1:]]
    ax_right.set_xticks(right_ticks)
    ax_right.set_xticklabels(right_labels, rotation=45, ha="right")
    ax_right.axvspan(0, n_inner, alpha=0.08)
    ax_right.axvline(n_inner, linestyle="--", linewidth=1)
    ax_right.set_xlim(0, RIGHT_TOTAL - 1)

    ax_left.spines["right"].set_visible(False)
    ax_mid.spines["left"].set_visible(False)
    ax_mid.spines["right"].set_visible(False)
    ax_right.spines["left"].set_visible(False)
    ax_mid.set_yticks([])
    ax_mid.tick_params(left=False)
    ax_right.yaxis.tick_right()
    ax_right.tick_params(labelright=False)

    kwargs = dict(color="k", clip_on=False, linewidth=1)
    ax_left.plot((1 - 0.015, 1 + 0.015), (-0.02, +0.02), transform=ax_left.transAxes, **kwargs)
    ax_left.plot((1 - 0.015, 1 + 0.015), (1 - 0.02, 1 + 0.02), transform=ax_left.transAxes, **kwargs)
    ax_mid.plot((-0.015, +0.015), (-0.02, +0.02), transform=ax_mid.transAxes, **kwargs)
    ax_mid.plot((-0.015, +0.015), (1 - 0.02, 1 + 0.02), transform=ax_mid.transAxes, **kwargs)
    ax_mid.plot((1 - 0.015, 1 + 0.015), (-0.02, +0.02), transform=ax_mid.transAxes, **kwargs)
    ax_mid.plot((1 - 0.015, 1 + 0.015), (1 - 0.02, 1 + 0.02), transform=ax_mid.transAxes, **kwargs)
    ax_right.plot((-0.015, +0.015), (-0.02, +0.02), transform=ax_right.transAxes, **kwargs)
    ax_right.plot((-0.015, +0.015), (1 - 0.02, 1 + 0.02), transform=ax_right.transAxes, **kwargs)


# -----------------------------------------------------
# plot: generate metagene
# -----------------------------------------------------
fig, (ax_left, ax_mid, ax_right) = plt.subplots(
    1, 3, figsize=(9, 4), sharey=True, gridspec_kw={"width_ratios": [LEFT_TOTAL, 2, RIGHT_TOTAL]}
)

left_x = np.arange(LEFT_TOTAL)
right_x = np.arange(RIGHT_TOTAL)
mid_x = np.array([0, 1])

label1 = Path(CO_FILE).name
ax_left.plot(left_x, co_left_raw, label=label1)
ax_mid.scatter(mid_x, [co_first_exon, co_last_exon], label=label1)
ax_right.plot(right_x, co_right_raw, label=label1)
if dsb_parts is not None:
    label2 = Path(DSB_FILE).name
    ax_left.plot(left_x, dsb_left_raw, label=label2)
    ax_mid.scatter(mid_x, [dsb_first_exon, dsb_last_exon], label=label2)
    ax_right.plot(right_x, dsb_right_raw, label=label2)
ax_left.set_ylabel("Density")
ax_right.legend()

style_dual_anchor(ax_left, ax_mid, ax_right)

plt.tight_layout()
plt.savefig("dual_anchor_metagene.pdf", bbox_inches="tight")
plt.show()
