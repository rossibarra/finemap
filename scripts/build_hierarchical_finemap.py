#!/usr/bin/env python3
"""Fit a smooth interval-censored recombination map.

For an observed crossover interval [L_i, R_i), the conditional likelihood is

    P(i | rate) = integral(L_i, R_i) rate(x) dx / integral(chrom) rate(x) dx.

The log-rate is constant within fixed-width bins.  A Gaussian random-walk
prior on adjacent log-rates supplies regularization.  The fitted shape is
finally scaled to the chromosome length of the Ogut genetic map.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).parent.parent
JRI = ROOT / "data/jri_v5.bed"
OGUT = ROOT / "data/ogut_fifthcM_map_agpv2.csv"
FAI = ROOT / "data/v5.fa.gz.fai"
DEFAULT_OUT = ROOT / "data/finemap_hierarchical_v5.bed"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=JRI)
    parser.add_argument("--ogut", type=Path, default=OGUT)
    parser.add_argument("--fai", type=Path, default=FAI)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument(
        "--smoothness",
        type=float,
        default=10.0,
        help="Gaussian random-walk precision on adjacent log-rates",
    )
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--learning-rate", type=float, default=0.03)
    parser.add_argument("--tolerance", type=float, default=1e-7)
    return parser.parse_args()


def interval_integrals(theta, edges, starts, ends):
    """Return rate integrals and bin indices for half-open intervals."""
    rate = np.exp(np.clip(theta, -30.0, 30.0))
    widths = np.diff(edges)
    cumulative = np.concatenate(([0.0], np.cumsum(rate * widths)))
    start_bin = np.minimum(starts // (edges[1] - edges[0]), len(rate) - 1)
    end_bin = np.minimum((ends - 1) // (edges[1] - edges[0]), len(rate) - 1)

    start_value = cumulative[start_bin] + rate[start_bin] * (starts - edges[start_bin])
    end_value = cumulative[end_bin] + rate[end_bin] * (ends - edges[end_bin])
    integrals = end_value - start_value
    return rate, widths, start_bin.astype(int), end_bin.astype(int), integrals


def expected_bin_counts(rate, edges, starts, ends, start_bin, end_bin, integrals):
    """Expected latent event allocations used by the likelihood gradient."""
    n_bins = len(rate)
    exposure_over_integral = np.zeros(n_bins)
    same = start_bin == end_bin

    if np.any(same):
        weights = (ends[same] - starts[same]) / integrals[same]
        exposure_over_integral += np.bincount(
            start_bin[same], weights=weights, minlength=n_bins
        )

    split = ~same
    if np.any(split):
        sb = start_bin[split]
        eb = end_bin[split]
        inv = 1.0 / integrals[split]
        left = (edges[sb + 1] - starts[split]) * inv
        right = (ends[split] - edges[eb]) * inv
        exposure_over_integral += np.bincount(sb, weights=left, minlength=n_bins)
        exposure_over_integral += np.bincount(eb, weights=right, minlength=n_bins)

        has_middle = eb > sb + 1
        difference = np.zeros(n_bins + 1)
        np.add.at(difference, sb[has_middle] + 1, inv[has_middle])
        np.add.at(difference, eb[has_middle], -inv[has_middle])
        exposure_over_integral += np.cumsum(difference[:-1]) * np.diff(edges)

    return rate * exposure_over_integral


def objective_and_gradient(theta, edges, starts, ends, smoothness):
    rate, widths, sb, eb, integrals = interval_integrals(theta, edges, starts, ends)
    if np.any(~np.isfinite(integrals)) or np.any(integrals <= 0):
        raise FloatingPointError("non-positive or non-finite interval likelihood")

    total = np.dot(rate, widths)
    n_events = len(starts)
    objective = np.log(integrals).sum() - n_events * np.log(total)
    allocated = expected_bin_counts(rate, edges, starts, ends, sb, eb, integrals)
    gradient = allocated - n_events * rate * widths / total

    differences = np.diff(theta)
    objective -= 0.5 * smoothness * np.dot(differences, differences)
    prior_gradient = np.zeros_like(theta)
    prior_gradient[:-1] += smoothness * differences
    prior_gradient[1:] -= smoothness * differences
    gradient += prior_gradient
    return objective, gradient


def fit_chromosome(starts, ends, chrom_length, bin_size, smoothness, iterations,
                   learning_rate, tolerance):
    edges = np.arange(0, chrom_length, bin_size, dtype=np.int64)
    edges = np.append(edges, chrom_length)
    theta = np.zeros(len(edges) - 1)
    first_moment = np.zeros_like(theta)
    second_moment = np.zeros_like(theta)
    previous = -np.inf

    for iteration in range(1, iterations + 1):
        objective, gradient = objective_and_gradient(
            theta, edges, starts, ends, smoothness
        )
        gradient /= len(starts)
        first_moment = 0.9 * first_moment + 0.1 * gradient
        second_moment = 0.999 * second_moment + 0.001 * gradient * gradient
        m_hat = first_moment / (1.0 - 0.9 ** iteration)
        v_hat = second_moment / (1.0 - 0.999 ** iteration)
        theta += learning_rate * m_hat / (np.sqrt(v_hat) + 1e-8)
        theta -= np.average(theta, weights=np.diff(edges))

        relative_change = (
            np.inf
            if not np.isfinite(previous)
            else abs(objective - previous) / max(1.0, abs(previous))
        )
        if iteration > 20 and relative_change < tolerance:
            break
        previous = objective

    final_objective, _ = objective_and_gradient(theta, edges, starts, ends, smoothness)
    return edges, np.exp(np.clip(theta, -30.0, 30.0)), iteration, final_objective


def load_inputs(args):
    jri = pd.read_csv(
        args.input, sep="\t", header=None, names=["chr", "start", "end", "id", "src"]
    )
    jri = jri[jri["end"] > jri["start"]].copy()
    ogut = pd.read_csv(args.ogut)
    ogut["chr"] = "Chr" + ogut["chromosome"].astype(str)
    targets = {
        chrom: group["cM"].max() - group["cM"].min()
        for chrom, group in ogut.groupby("chr")
    }
    fai = pd.read_csv(args.fai, sep="\t", header=None, usecols=[0, 1], names=["seq", "length"])
    fai["chr"] = fai["seq"].str.replace(r"^chr", "Chr", regex=True)
    lengths = dict(zip(fai["chr"], fai["length"]))
    return jri, targets, lengths


def main():
    args = parse_args()
    if args.bin_size <= 0 or args.smoothness < 0 or args.iterations <= 0:
        raise SystemExit("bin size and iterations must be positive; smoothness must be nonnegative")

    jri, targets, lengths = load_inputs(args)
    records = []
    for chrom in sorted(jri["chr"].unique(), key=lambda value: int(value[3:])):
        if chrom not in targets or chrom not in lengths:
            print(f"Skipping {chrom}: missing length or Ogut target", file=sys.stderr)
            continue
        subset = jri[jri["chr"] == chrom]
        starts = subset["start"].to_numpy(dtype=np.int64)
        ends = subset["end"].to_numpy(dtype=np.int64)
        chrom_length = int(lengths[chrom])
        if np.any(starts < 0) or np.any(ends > chrom_length):
            raise SystemExit(f"interval outside chromosome bounds on {chrom}")

        edges, relative_rate, iterations, objective = fit_chromosome(
            starts, ends, chrom_length, args.bin_size, args.smoothness,
            args.iterations, args.learning_rate, args.tolerance,
        )
        widths = np.diff(edges)
        cM_per_bp = relative_rate * targets[chrom] / np.dot(relative_rate, widths)
        cumulative = np.concatenate(([0.0], np.cumsum(cM_per_bp * widths)))
        for index, rate in enumerate(cM_per_bp):
            records.append((
                chrom, int(edges[index]), int(edges[index + 1]),
                cumulative[index], cumulative[index + 1], rate * 1e6,
            ))
        print(
            f"{chrom}: {len(starts)} intervals, {len(relative_rate)} bins, "
            f"{iterations} iterations, objective={objective:.6f}",
            file=sys.stderr,
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        for record in records:
            handle.write("\t".join(map(str, record)) + "\n")
    print(f"Wrote {len(records)} segments to {args.output}")


if __name__ == "__main__":
    main()
