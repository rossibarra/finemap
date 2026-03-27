#!/usr/bin/env python3

#### module load python/3.9.19

#### USAGE:
####   script.py --gff genes.gff3 --bed1 co.bed --bed2 dsb.bed --flank 10000 20 --body 50 --gene 2000 20
####   script.py --gene_bed genes.bed --bed1 co.bed --flank 10000 20 --body 50 --gene 2000 20

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
    parser.add_argument("--body", type=int, required=True, help="Number of bins for the scaled gene body.")
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
BODY_BIN = args.body
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
# define gene 
# =====================================================

n_up    = FLANK_N_BINS
n_inner = INNER_N_BINS
n_body  = BODY_BIN  

TOTAL = n_up + n_inner + n_body + n_inner + n_up

# offsets
OFF_TSS  = n_up
OFF_BODY = OFF_TSS + n_inner
OFF_TTS  = OFF_BODY + n_body
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

    profile = np.zeros(TOTAL)
    if g.chr not in signal_dict:
        return profile

    if UNIFORM:
        starts, ends, weights = signal_dict[g.chr]

        def range_signal(lo, hi):
            left = max(0, np.searchsorted(starts, lo, side="right") - 1)
            right = np.searchsorted(starts, hi, side="left")
            total = 0.0
            for idx in range(left, right):
                ov_start = max(lo, starts[idx])
                ov_end = min(hi, ends[idx])
                if ov_end > ov_start:
                    total += (ov_end - ov_start) * weights[idx]
            return total
    else:
        m, c = signal_dict[g.chr]

        def slice_range(lo, hi):
            left = np.searchsorted(m, lo, side="left")
            right = np.searchsorted(m, hi, side="left")
            return m[left:right], c[left:right]

    # strand normalization
    if g.strand == "+":
        TSS, TTS = g.start, g.end
    else:
        TSS, TTS = g.end, g.start

    # gene body limits
    body_start = TSS + INNER_SIZE
    body_end   = TTS - INNER_SIZE
    body_len   = max(1, body_end - body_start)

    # --------------- upstream ---------------
    if UNIFORM:
        for i in range(n_up):
            b0 = TSS - FLANK_SIZE + i * FLANK_SIZE / n_up
            b1 = TSS - FLANK_SIZE + (i + 1) * FLANK_SIZE / n_up
            profile[i] += range_signal(int(round(b0)), int(round(b1)))
    else:
        m_sub, c_sub = slice_range(TSS - FLANK_SIZE, TSS)
        idx = (((m_sub - (TSS - FLANK_SIZE)) / FLANK_SIZE) * n_up).astype(int)
        valid = (idx >= 0) & (idx < n_up)
        np.add.at(profile, idx[valid], c_sub[valid])

    # --------------- TSS inner ---------------
    if UNIFORM:
        for i in range(n_inner):
            b0 = TSS + i * INNER_SIZE / n_inner
            b1 = TSS + (i + 1) * INNER_SIZE / n_inner
            profile[OFF_TSS + i] += range_signal(int(round(b0)), int(round(b1)))
    else:
        m_sub, c_sub = slice_range(TSS, TSS + INNER_SIZE)
        idx = OFF_TSS + (((m_sub - TSS) / INNER_SIZE) * n_inner).astype(int)
        valid = (idx >= OFF_TSS) & (idx < OFF_BODY)
        np.add.at(profile, idx[valid], c_sub[valid])

    # --------------- scaled gene body ---------------
    if UNIFORM:
        for i in range(n_body):
            b0 = body_start + i * body_len / n_body
            b1 = body_start + (i + 1) * body_len / n_body
            profile[OFF_BODY + i] += range_signal(int(round(b0)), int(round(b1)))
    else:
        m_sub, c_sub = slice_range(body_start, body_end)
        scaled = ((m_sub - body_start) / body_len) * n_body
        idx = OFF_BODY + scaled.astype(int)
        valid = (idx >= OFF_BODY) & (idx < OFF_TTS)
        np.add.at(profile, idx[valid], c_sub[valid])

    # --------------- TTS inner ---------------
    if UNIFORM:
        for i in range(n_inner):
            b0 = TTS - INNER_SIZE + i * INNER_SIZE / n_inner
            b1 = TTS - INNER_SIZE + (i + 1) * INNER_SIZE / n_inner
            profile[OFF_TTS + i] += range_signal(int(round(b0)), int(round(b1)))
    else:
        m_sub, c_sub = slice_range(TTS - INNER_SIZE, TTS)
        idx = OFF_TTS + (((m_sub - (TTS - INNER_SIZE)) / INNER_SIZE) * n_inner).astype(int)
        valid = (idx >= OFF_TTS) & (idx < OFF_DOWN)
        np.add.at(profile, idx[valid], c_sub[valid])

    # --------------- downstream ---------------
    if UNIFORM:
        for i in range(n_up):
            b0 = TTS + i * FLANK_SIZE / n_up
            b1 = TTS + (i + 1) * FLANK_SIZE / n_up
            profile[OFF_DOWN + i] += range_signal(int(round(b0)), int(round(b1)))
    else:
        m_sub, c_sub = slice_range(TTS, TTS + FLANK_SIZE)
        idx = OFF_DOWN + (((m_sub - TTS) / FLANK_SIZE) * n_up).astype(int)
        valid = (idx >= OFF_DOWN) & (idx < TOTAL)
        np.add.at(profile, idx[valid], c_sub[valid])

    return profile

# =====================================================
# aggregation function
# =====================================================

def aggregate(signal_dict):

    profiles = []

    for g in genes.itertuples(index=False):
        p = gene_profile(g, signal_dict)
        if p.sum() == 0:
            continue
        profiles.append(p)

    if not profiles:
        raise SystemExit(
            "No overlapping signal was found between the gene annotations and signal BED files after chromosome normalization."
        )
    profiles = np.vstack(profiles)
    return profiles.mean(axis=0) #metagene profile

# =====================================================
# smoothing function
# =====================================================

def smooth(x, w=3):
    return np.convolve(x, np.ones(w)/w, mode="same")

# =====================================================
# run aggregation & smoothing
# =====================================================

co_raw = aggregate(co_chr)
dsb_raw = aggregate(dsb_chr) if dsb_chr is not None else None

co_raw  = smooth(co_raw, SMOOTH)
if dsb_raw is not None:
    dsb_raw = smooth(dsb_raw, SMOOTH)

# -----------------------------------------------------
# plot: define boundaries
# -----------------------------------------------------

TSS_START   = OFF_TSS
TSS_END     = OFF_BODY

BODY_START  = OFF_BODY
BODY_END    = OFF_TTS

TTS_START   = OFF_TTS
TTS_END     = OFF_DOWN

BODY_CENTER = (BODY_START + BODY_END) / 2

# -----------------------------------------------------
# plot: helper function to format x-axis
# -----------------------------------------------------

def style_dual_anchor(ax):

    # Shade TSS and TTS windows
    ax.axvspan(TSS_START, TSS_END, alpha=0.08)
    ax.axvspan(TTS_START, TTS_END, alpha=0.08)

    # Anchor TSS and TTS
    ax.axvline(TSS_START, linestyle="--", linewidth=1)
    ax.axvline(TTS_END, linestyle="--", linewidth=1)

    # Clean x-axis labels
    xticks = [
        0,
        TSS_START,
        BODY_CENTER,
        TTS_END,
        TOTAL - 1
    ]

    xlabels = [
        "-10 kb",
        "TSS",
        "Gene body",
        "TTS",
        "+10 kb"
    ]

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_xlim(0, TOTAL - 1)


# -----------------------------------------------------
# plot: generate metagene
# -----------------------------------------------------
fig, ax = plt.subplots(1, 1, figsize=(6, 4))

x = np.arange(TOTAL)

label1 = Path(CO_FILE).name
ax.plot(x, co_raw, label=label1)
if dsb_raw is not None:
    ax.plot(x, dsb_raw, label=Path(DSB_FILE).name)
ax.set_ylabel("Density")
ax.legend()

style_dual_anchor(ax)

plt.tight_layout()
plt.savefig("dual_anchor_metagene.pdf", bbox_inches="tight")
plt.show()
