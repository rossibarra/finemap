#!/usr/bin/env python3

#### module load python/3.9.19

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot dual-anchor metagene profiles for one or two signal BED files."
    )
    gene_group = parser.add_mutually_exclusive_group(required=True)
    gene_group.add_argument(
        "--gff3",
        dest="gff3_file",
        help="Gene annotation in GFF/GFF3 format.",
    )
    gene_group.add_argument(
        "--gene-bed",
        dest="gene_bed_file",
        help="Legacy BED-like gene file with columns: chr, start, end, len, dot, strand, type.",
    )
    parser.add_argument(
        "--co-bed",
        required=True,
        help="CO signal BED with columns: chr, start, end, count.",
    )
    parser.add_argument(
        "--dsb-bed",
        help="Optional DSB signal BED with columns: chr, start, end, count.",
    )
    parser.add_argument(
        "--flank-bin-size",
        type=int,
        required=True,
        help="Bin size in bp for upstream and downstream flanks.",
    )
    parser.add_argument(
        "--body-bins",
        type=int,
        required=True,
        help="Number of bins used for the scaled gene body.",
    )
    parser.add_argument(
        "--inner-bin-size",
        type=int,
        required=True,
        help="Bin size in bp for the TSS and TTS inner windows.",
    )
    parser.add_argument(
        "--flank-size",
        type=int,
        default=10000,
        help="Flank size in bp on each side of the gene. Default: 10000.",
    )
    parser.add_argument(
        "--inner-size",
        type=int,
        default=2000,
        help="Inner window size in bp at TSS and TTS. Default: 2000.",
    )
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=5,
        help="Moving-average smoothing window size. Default: 5.",
    )
    parser.add_argument(
        "--output",
        default="dual_anchor_metagene.pdf",
        help="Output PDF path. Default: dual_anchor_metagene.pdf.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Write the PDF without opening an interactive plot window.",
    )
    return parser.parse_args()


ARGS = parse_args()

GENE_FILE = ARGS.gff3_file or ARGS.gene_bed_file
CO_FILE = ARGS.co_bed
DSB_FILE = ARGS.dsb_bed
FLANK_BIN = ARGS.flank_bin_size
BODY_BIN = ARGS.body_bins
INNER_BIN = ARGS.inner_bin_size
FLANK_SIZE = ARGS.flank_size
INNER_SIZE = ARGS.inner_size
SMOOTH = ARGS.smooth_window
OUTPUT = ARGS.output


def normalize_chrom_name(value):
    if pd.isna(value):
        return value
    text = str(value).strip()
    if text.lower().startswith("chr"):
        return f"chr{text[3:]}"
    return text


def load_genes(path):
    if path.endswith((".gff3", ".gff", ".gff3.gz", ".gff.gz")):
        genes = pd.read_csv(
            path,
            sep="\t",
            header=None,
            comment="#",
            names=[
                "chr",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ],
            usecols=list(range(9)),
        )
        genes = genes[genes["feature"] == "gene"].copy()
        genes = genes[["chr", "start", "end", "strand"]]
        genes["start"] = pd.to_numeric(genes["start"], errors="coerce") - 1
        genes["end"] = pd.to_numeric(genes["end"], errors="coerce")
        genes["len"] = genes["end"] - genes["start"]
        genes["dot"] = "."
        genes["type"] = "gene"
        genes = genes[["chr", "start", "end", "len", "dot", "strand", "type"]]
    else:
        genes = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "len", "dot", "strand", "type"],
            usecols=list(range(7)),
        )

    genes["chr"] = genes["chr"].map(normalize_chrom_name)
    genes["start"] = pd.to_numeric(genes["start"], errors="coerce")
    genes["end"] = pd.to_numeric(genes["end"], errors="coerce")
    genes = genes.dropna(subset=["chr", "start", "end", "strand"]).copy()
    genes["start"] = genes["start"].astype(np.int64)
    genes["end"] = genes["end"].astype(np.int64)
    genes["strand"] = genes["strand"].astype(str)
    return genes


def load_signal(path):
    signal = pd.read_csv(path, sep="\t", header=None)
    if signal.shape[1] < 3:
        raise SystemExit(f"Signal BED must have at least 3 columns: {path}")

    signal = signal.iloc[:, :4].copy()
    signal.columns = ["chr", "start", "end", "count"][: signal.shape[1]]
    if "count" not in signal.columns:
        signal["count"] = 1.0

    signal["chr"] = signal["chr"].map(normalize_chrom_name)
    signal["start"] = pd.to_numeric(signal["start"], errors="coerce")
    signal["end"] = pd.to_numeric(signal["end"], errors="coerce")
    signal["count"] = pd.to_numeric(signal["count"], errors="coerce")
    signal = signal.dropna(subset=["chr", "start", "end"]).copy()
    signal["count"] = signal["count"].fillna(1.0)
    signal["start"] = signal["start"].astype(np.int64)
    signal["end"] = signal["end"].astype(np.int64)
    signal["mid"] = ((signal["start"] + signal["end"]) // 2).astype(np.int64)
    return signal


genes = load_genes(GENE_FILE)
co = load_signal(CO_FILE)
dsb = load_signal(DSB_FILE) if DSB_FILE else None

n_up = FLANK_SIZE // FLANK_BIN
n_inner = INNER_SIZE // INNER_BIN
n_body = BODY_BIN

TOTAL = n_up + n_inner + n_body + n_inner + n_up

OFF_TSS = n_up
OFF_BODY = OFF_TSS + n_inner
OFF_TTS = OFF_BODY + n_body
OFF_DOWN = OFF_TTS + n_inner


def build_gene_dict(df):
    gene_dict = {}
    for chrom, sub in df.groupby("chr", sort=False):
        gene_dict[chrom] = {
            "start": sub["start"].to_numpy(dtype=np.int64, copy=True),
            "end": sub["end"].to_numpy(dtype=np.int64, copy=True),
            "strand": sub["strand"].to_numpy(copy=True),
        }
    return gene_dict


def build_signal_dict(df):
    signal_dict = {}
    for chrom, sub in df.groupby("chr", sort=False):
        mids = sub["mid"].to_numpy(dtype=np.int64, copy=True)
        counts = sub["count"].to_numpy(dtype=np.float64, copy=True)
        order = np.argsort(mids, kind="mergesort")
        signal_dict[chrom] = {
            "mid": mids[order],
            "count": counts[order],
        }
    return signal_dict


GENE_DICT = build_gene_dict(genes)
CO_DICT = build_signal_dict(co)
DSB_DICT = build_signal_dict(dsb) if dsb is not None else None


def add_binned_counts(profile, mids, counts, low, high, offset, bin_size, lower_bound, upper_bound):
    left = np.searchsorted(mids, low, side="left")
    right = np.searchsorted(mids, high, side="left")
    region_mids = mids[left:right]
    if region_mids.size == 0:
        return
    idx = offset + ((region_mids - low) // bin_size).astype(np.int64)
    valid = (idx >= lower_bound) & (idx < upper_bound)
    if np.any(valid):
        # Normalize fixed-width regions to density per bp.
        np.add.at(profile, idx[valid], counts[left:right][valid] / float(bin_size))


def add_scaled_body_counts(profile, mids, counts, body_start, body_end):
    left = np.searchsorted(mids, body_start, side="left")
    right = np.searchsorted(mids, body_end, side="left")
    region_mids = mids[left:right]
    if region_mids.size == 0:
        return
    body_len = max(1, body_end - body_start)
    scaled = ((region_mids - body_start) / body_len) * n_body
    idx = OFF_BODY + scaled.astype(np.int64)
    valid = (idx >= OFF_BODY) & (idx < OFF_TTS)
    if np.any(valid):
        # Gene-body bins have gene-specific widths, so normalize by that width.
        body_bin_width = body_len / float(n_body)
        np.add.at(profile, idx[valid], counts[left:right][valid] / body_bin_width)


def gene_profile(chrom, start, end, strand, signal_dict):
    profile = np.zeros(TOTAL, dtype=np.float64)
    chrom_signal = signal_dict.get(chrom)
    if chrom_signal is None:
        return profile

    mids = chrom_signal["mid"]
    counts = chrom_signal["count"]

    if strand == "+":
        tss = start
        tts = end
    else:
        tss = end
        tts = start

    body_start = tss + INNER_SIZE
    body_end = tts - INNER_SIZE

    add_binned_counts(
        profile,
        mids,
        counts,
        tss - FLANK_SIZE,
        tss,
        0,
        FLANK_BIN,
        0,
        n_up,
    )
    add_binned_counts(
        profile,
        mids,
        counts,
        tss,
        tss + INNER_SIZE,
        OFF_TSS,
        INNER_BIN,
        OFF_TSS,
        OFF_BODY,
    )
    add_scaled_body_counts(profile, mids, counts, body_start, body_end)
    add_binned_counts(
        profile,
        mids,
        counts,
        tts - INNER_SIZE,
        tts,
        OFF_TTS,
        INNER_BIN,
        OFF_TTS,
        OFF_DOWN,
    )
    add_binned_counts(
        profile,
        mids,
        counts,
        tts,
        tts + FLANK_SIZE,
        OFF_DOWN,
        FLANK_BIN,
        OFF_DOWN,
        TOTAL,
    )
    return profile


def aggregate(signal_dict):
    profile_sum = np.zeros(TOTAL, dtype=np.float64)
    n_profiles = 0

    for chrom, gene_arrays in GENE_DICT.items():
        if chrom not in signal_dict:
            continue
        starts = gene_arrays["start"]
        ends = gene_arrays["end"]
        strands = gene_arrays["strand"]
        for start, end, strand in zip(starts, ends, strands):
            profile = gene_profile(chrom, start, end, strand, signal_dict)
            if np.any(profile):
                profile_sum += profile
                n_profiles += 1

    if n_profiles == 0:
        return profile_sum
    return profile_sum / n_profiles


def smooth(values, window):
    if window <= 1:
        return values
    kernel = np.ones(window, dtype=np.float64) / window
    return np.convolve(values, kernel, mode="same")


co_raw = smooth(aggregate(CO_DICT), SMOOTH)
dsb_raw = smooth(aggregate(DSB_DICT), SMOOTH) if DSB_DICT is not None else None

TSS_START = OFF_TSS
TSS_END = OFF_BODY
BODY_START = OFF_BODY
BODY_END = OFF_TTS
TTS_START = OFF_TTS
TTS_END = OFF_DOWN
BODY_CENTER = (BODY_START + BODY_END) / 2


def style_dual_anchor(ax):
    ax.axvspan(TSS_START, TSS_END, alpha=0.08)
    ax.axvspan(TTS_START, TTS_END, alpha=0.08)
    ax.axvline(TSS_START, linestyle="--", linewidth=1)
    ax.axvline(TTS_END, linestyle="--", linewidth=1)
    ax.set_xticks([0, TSS_START, BODY_CENTER, TTS_END, TOTAL - 1])
    flank_kb = FLANK_SIZE / 1000.0
    if flank_kb.is_integer():
        flank_label = f"{int(flank_kb)} kb"
    else:
        flank_label = f"{flank_kb:g} kb"
    ax.set_xticklabels([f"-{flank_label}", "TSS", "Gene body", "TTS", f"+{flank_label}"])
    ax.set_xlim(0, TOTAL - 1)


fig, ax = plt.subplots(1, 1, figsize=(6, 4))
x = np.arange(TOTAL)

ax.plot(x, co_raw, label="CO")
if dsb_raw is not None:
    ax.plot(x, dsb_raw, label="DSB")
ax.set_title("Dual-anchor metagene")
ax.set_ylabel("Density per bp")
ax.legend()

style_dual_anchor(ax)

plt.tight_layout()
plt.savefig(OUTPUT, bbox_inches="tight")
if not ARGS.no_show:
    plt.show()
else:
    plt.close(fig)
