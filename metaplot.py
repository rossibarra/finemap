#!/usr/bin/env python3

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


WARNED_CHROM_RENAMES = set()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a TSS/TTS metaplot from GFF gene annotations and a signal BED."
    )
    parser.add_argument("--gff", required=True, help="Gene annotation GFF/GFF3 file.")
    parser.add_argument(
        "--input",
        required=True,
        help="Signal BED file. Midpoints are assigned to bins; missing or non-numeric column 4 is treated as 1.",
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        required=True,
        help="Bin size in bp.",
    )
    parser.add_argument(
        "--flanking-bp",
        type=int,
        required=True,
        help="Maximum flank size in bp to include upstream of TSS and downstream of TTS.",
    )
    parser.add_argument(
        "--body-bins",
        type=int,
        required=True,
        help="Number of internal bins to plot on each side of the gene body.",
    )
    parser.add_argument(
        "--output",
        default="metaplot.pdf",
        help="Output PDF path. Default: metaplot.pdf.",
    )
    parser.add_argument(
        "--title",
        default="Metaplot",
        help="Plot title. Default: Metaplot.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Write the plot without opening an interactive window.",
    )
    return parser.parse_args()


ARGS = parse_args()

if ARGS.bin_size <= 0:
    raise SystemExit("--bin-size must be > 0")
if ARGS.flanking_bp < 0:
    raise SystemExit("--flanking-bp must be >= 0")
if ARGS.flanking_bp % ARGS.bin_size != 0:
    raise SystemExit("--flanking-bp must be divisible by --bin-size")
if ARGS.body_bins <= 0:
    raise SystemExit("--body-bins must be > 0")


def normalize_chrom_name(value):
    if pd.isna(value):
        return value
    text = str(value).strip()
    if text.lower().startswith("chr"):
        normalized = f"chr{text[3:]}"
    else:
        normalized = text
    if normalized != text and text not in WARNED_CHROM_RENAMES:
        print(
            f"Warning: normalized chromosome name '{text}' -> '{normalized}'",
            file=sys.stderr,
        )
        WARNED_CHROM_RENAMES.add(text)
    return normalized


def load_genes(path):
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
    genes["chr"] = genes["chr"].map(normalize_chrom_name)
    genes["start"] = pd.to_numeric(genes["start"], errors="coerce") - 1
    genes["end"] = pd.to_numeric(genes["end"], errors="coerce")
    genes = genes.dropna(subset=["chr", "start", "end", "strand"]).copy()
    genes["start"] = genes["start"].astype(np.int64)
    genes["end"] = genes["end"].astype(np.int64)
    genes["strand"] = genes["strand"].astype(str)
    genes = genes[genes["strand"].isin(["+", "-"])].copy()
    return genes


def load_signal(path):
    signal = pd.read_csv(path, sep="\t", header=None)
    if signal.shape[1] < 3:
        raise SystemExit(f"Signal BED must have at least 3 columns: {path}")
    if signal.shape[1] >= 4:
        signal = signal.iloc[:, :4].copy()
        signal.columns = ["chr", "start", "end", "value"]
    else:
        signal = signal.iloc[:, :3].copy()
        signal.columns = ["chr", "start", "end"]
        signal["value"] = 1.0

    signal["chr"] = signal["chr"].map(normalize_chrom_name)
    signal["start"] = pd.to_numeric(signal["start"], errors="coerce")
    signal["end"] = pd.to_numeric(signal["end"], errors="coerce")
    signal["value"] = pd.to_numeric(signal["value"], errors="coerce").fillna(1.0)
    signal = signal.dropna(subset=["chr", "start", "end"]).copy()
    signal["start"] = signal["start"].astype(np.int64)
    signal["end"] = signal["end"].astype(np.int64)
    signal["mid"] = ((signal["start"] + signal["end"]) // 2).astype(np.int64)
    return signal


def build_signal_dict(df):
    signal_dict = {}
    for chrom, sub in df.groupby("chr", sort=False):
        mids = sub["mid"].to_numpy(dtype=np.int64, copy=True)
        values = sub["value"].to_numpy(dtype=np.float64, copy=True)
        order = np.argsort(mids, kind="mergesort")
        signal_dict[chrom] = {
            "mid": mids[order],
            "value": values[order],
        }
    return signal_dict


def midpoint_sum(chrom_signal, window_start, window_end):
    mids = chrom_signal["mid"]
    values = chrom_signal["value"]
    left = np.searchsorted(mids, window_start, side="left")
    right = np.searchsorted(mids, window_end, side="left")
    if right <= left:
        return 0.0
    return float(values[left:right].sum())


def add_window(signal_dict, chrom, slot, window_start, window_end, sums, sumsq, counts):
    chrom_signal = signal_dict.get(chrom)
    if chrom_signal is None:
        value = 0.0
    else:
        value = midpoint_sum(chrom_signal, window_start, window_end)
    sums[slot] += value
    sumsq[slot] += value * value
    counts[slot] += 1


def window_stats(sums, sumsq, counts):
    means = np.divide(
        sums,
        counts,
        out=np.full_like(sums, np.nan),
        where=counts > 0,
    )
    variances = np.divide(
        sumsq,
        counts,
        out=np.zeros_like(sums),
        where=counts > 0,
    ) - np.square(np.nan_to_num(means, nan=0.0))
    variances = np.maximum(variances, 0.0)
    ses = np.full_like(sums, np.nan)
    valid = counts > 0
    ses[valid] = np.sqrt(variances[valid] / counts[valid])
    return means, ses


def iter_plus_windows(start, end, bin_size):
    pos = start
    while pos + bin_size <= end:
        yield int(pos), int(pos + bin_size)
        pos += bin_size


def iter_minus_windows(start, end, bin_size):
    pos = end
    while pos - bin_size >= start:
        yield int(pos - bin_size), int(pos)
        pos -= bin_size


def prepare_genes(genes):
    prepared = []
    for chrom, sub in genes.groupby("chr", sort=False):
        sub = sub.sort_values(["start", "end"], kind="mergesort").reset_index(drop=True)
        starts = sub["start"].to_numpy(dtype=np.int64, copy=True)
        ends = sub["end"].to_numpy(dtype=np.int64, copy=True)
        prev_end = np.full(len(sub), -1, dtype=np.int64)
        next_start = np.full(len(sub), np.iinfo(np.int64).max, dtype=np.int64)
        if len(sub) > 1:
            prev_end[1:] = np.maximum.accumulate(ends[:-1])
            next_start[:-1] = starts[1:]

        for idx, row in sub.iterrows():
            strand = row["strand"]
            start = int(row["start"])
            end = int(row["end"])
            if strand == "+":
                tss = start
                tts = end
                flank5_low = max(tss - ARGS.flanking_bp, int(prev_end[idx]))
                flank5_high = tss
                flank3_low = tts
                flank3_high = min(tts + ARGS.flanking_bp, int(next_start[idx]))
                gene_oriented_start = start
                gene_oriented_end = end
                iter_gene_from_tss = iter_plus_windows
                iter_gene_to_tts = iter_minus_windows
                flank5_reverse = False
                flank3_reverse = False
            else:
                tss = end
                tts = start
                flank5_low = tss
                flank5_high = min(tss + ARGS.flanking_bp, int(next_start[idx]))
                flank3_low = max(tts - ARGS.flanking_bp, int(prev_end[idx]))
                flank3_high = tts
                gene_oriented_start = start
                gene_oriented_end = end
                iter_gene_from_tss = iter_minus_windows
                iter_gene_to_tts = iter_plus_windows
                flank5_reverse = True
                flank3_reverse = True

            prepared.append(
                {
                    "chr": chrom,
                    "strand": strand,
                    "tss": tss,
                    "tts": tts,
                    "gene_start": gene_oriented_start,
                    "gene_end": gene_oriented_end,
                    "flank5_low": int(flank5_low),
                    "flank5_high": int(flank5_high),
                    "flank3_low": int(flank3_low),
                    "flank3_high": int(flank3_high),
                    "iter_gene_from_tss": iter_gene_from_tss,
                    "iter_gene_to_tts": iter_gene_to_tts,
                    "flank5_reverse": flank5_reverse,
                    "flank3_reverse": flank3_reverse,
                }
            )
    return prepared


def aggregate_profiles(genes, signal_dict):
    flank_bins = ARGS.flanking_bp // ARGS.bin_size
    total_bins = flank_bins + (2 * ARGS.body_bins) + flank_bins
    sums = np.zeros(total_bins, dtype=np.float64)
    sumsq = np.zeros(total_bins, dtype=np.float64)
    counts = np.zeros(total_bins, dtype=np.int64)

    left_internal_offset = flank_bins
    right_internal_offset = flank_bins + ARGS.body_bins
    right_flank_offset = flank_bins + (2 * ARGS.body_bins)

    for gene in genes:
        chrom = gene["chr"]

        flank5_bins = list(iter_plus_windows(gene["flank5_low"], gene["flank5_high"], ARGS.bin_size))
        if gene["flank5_reverse"]:
            flank5_bins.reverse()
        flank5_slots_start = flank_bins - len(flank5_bins)
        for local_idx, (window_start, window_end) in enumerate(flank5_bins):
            slot = flank5_slots_start + local_idx
            add_window(signal_dict, chrom, slot, window_start, window_end, sums, sumsq, counts)

        full_gene_bins = list(gene["iter_gene_from_tss"](gene["gene_start"], gene["gene_end"], ARGS.bin_size))
        total_gene_bins = len(full_gene_bins)
        if total_gene_bins <= 2 * ARGS.body_bins:
            per_side_bins = total_gene_bins // 2
            left_gene_bins = full_gene_bins[:per_side_bins]
            right_gene_bins = list(
                gene["iter_gene_to_tts"](gene["gene_start"], gene["gene_end"], ARGS.bin_size)
            )[:per_side_bins]
        else:
            left_gene_bins = full_gene_bins[:ARGS.body_bins]
            right_gene_bins = list(
                gene["iter_gene_to_tts"](gene["gene_start"], gene["gene_end"], ARGS.bin_size)
            )[:ARGS.body_bins]
        right_gene_bins.reverse()

        for local_idx, (window_start, window_end) in enumerate(left_gene_bins):
            slot = left_internal_offset + local_idx
            add_window(signal_dict, chrom, slot, window_start, window_end, sums, sumsq, counts)

        for local_idx, (window_start, window_end) in enumerate(right_gene_bins):
            slot = right_internal_offset + local_idx
            add_window(signal_dict, chrom, slot, window_start, window_end, sums, sumsq, counts)

        flank3_bins = list(iter_plus_windows(gene["flank3_low"], gene["flank3_high"], ARGS.bin_size))
        if gene["flank3_reverse"]:
            flank3_bins.reverse()
        for local_idx, (window_start, window_end) in enumerate(flank3_bins):
            slot = right_flank_offset + local_idx
            add_window(signal_dict, chrom, slot, window_start, window_end, sums, sumsq, counts)

    averages, ses = window_stats(sums, sumsq, counts)
    return averages, ses, counts


def build_plot(averages, ses):
    flank_bins = ARGS.flanking_bp // ARGS.bin_size
    total_bins = len(averages)
    x = np.arange(total_bins, dtype=np.int64)

    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    ax.fill_between(x, averages - ses, averages + ses, alpha=0.25, linewidth=0)
    ax.plot(x, averages, linewidth=1.5)
    ax.axvline(flank_bins - 0.5, linestyle="--", linewidth=1)
    ax.axvline(flank_bins + (2 * ARGS.body_bins) - 0.5, linestyle="--", linewidth=1)
    ax.set_ylabel("Average signal")
    ax.set_title(ARGS.title)

    internal_bp = 2 * ARGS.body_bins * ARGS.bin_size
    flank_kb = ARGS.flanking_bp / 1000.0
    internal_kb = internal_bp / 1000.0

    def kb_label(value):
        if float(value).is_integer():
            return f"{int(value)} kb"
        return f"{value:g} kb"

    ax.set_xticks([0, flank_bins, flank_bins + (2 * ARGS.body_bins), total_bins - 1])
    ax.set_xticklabels(
        [
            f"-{kb_label(flank_kb)}",
            "TSS",
            "TTS",
            f"+{kb_label(flank_kb)}",
        ]
    )
    ax.set_xlabel(f"Flank / gene window ({kb_label(internal_kb)} shown inside gene)")
    ax.set_xlim(0, total_bins - 1)
    plt.tight_layout()
    return fig


GENES = load_genes(ARGS.gff)
SIGNAL = load_signal(ARGS.input)
GENE_LAYOUT = prepare_genes(GENES)
SIGNAL_DICT = build_signal_dict(SIGNAL)
AVERAGES, SES, COUNTS = aggregate_profiles(GENE_LAYOUT, SIGNAL_DICT)

FIG = build_plot(AVERAGES, SES)
FIG.savefig(ARGS.output, bbox_inches="tight")
if not ARGS.no_show:
    plt.show()
else:
    plt.close(FIG)
