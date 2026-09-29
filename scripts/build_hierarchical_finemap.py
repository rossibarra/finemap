#!/usr/bin/env python3
"""Fit a smooth interval-censored recombination map.

For an observed crossover interval [L_i, R_i), the conditional likelihood is

    P(i | rate) = integral(L_i, R_i) rate(x) dx / integral(chrom) rate(x) dx.

The log-rate is constant within fixed-width bins.  A Gaussian random-walk
prior on adjacent log-rates supplies regularization.  The fitted shape is
finally scaled to the chromosome length of the Ogut genetic map.

The penalized likelihood is maximized per chromosome by a line-searched
Newton-CG method that uses exact Hessian-vector products.  A fit is accepted
only when both (i) the largest absolute gradient of the log posterior with
respect to any bin's log-rate is below --tolerance and (ii) no normalized bin
rate has changed by more than --rate-tolerance (relative) over the last
--stability-window iterations.  Reaching --iterations first is an error.
"""

import argparse
import sys
import time
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
    parser.add_argument("--iterations", type=int, default=200,
                        help="maximum Newton iterations per chromosome; hitting it is an error")
    parser.add_argument("--learning-rate", type=float, default=1.0,
                        help="initial Newton line-search step (backtracked as needed)")
    parser.add_argument("--tolerance", type=float, default=1e-6,
                        help="convergence threshold on max |d log posterior / d log-rate| "
                             "over bins (units: expected crossovers per bin)")
    parser.add_argument("--rate-tolerance", type=float, default=1e-6,
                        help="convergence threshold on the max relative change of any "
                             "normalized bin rate over the last --stability-window iterations")
    parser.add_argument("--stability-window", type=int, default=3,
                        help="number of iterations over which rate stability is assessed")
    return parser.parse_args()


def interval_integrals(theta, edges, starts, ends):
    """Return rate integrals and bin indices for half-open intervals."""
    rate = np.exp(theta)
    widths = np.diff(edges)
    cumulative = np.concatenate(([0.0], np.cumsum(rate * widths)))
    start_bin = np.minimum(starts // (edges[1] - edges[0]), len(rate) - 1)
    end_bin = np.minimum((ends - 1) // (edges[1] - edges[0]), len(rate) - 1)

    start_value = cumulative[start_bin] + rate[start_bin] * (starts - edges[start_bin])
    end_value = cumulative[end_bin] + rate[end_bin] * (ends - edges[end_bin])
    integrals = end_value - start_value
    return rate, widths, start_bin.astype(int), end_bin.astype(int), integrals


def back_project(rate, edges, starts, ends, start_bin, end_bin, coefficients):
    """Return rate_j * sum_i overlap(i, j) * coefficients_i for every bin j."""
    n_bins = len(rate)
    accumulated = np.zeros(n_bins)
    same = start_bin == end_bin

    if np.any(same):
        weights = (ends[same] - starts[same]) * coefficients[same]
        accumulated += np.bincount(start_bin[same], weights=weights, minlength=n_bins)

    split = ~same
    if np.any(split):
        sb = start_bin[split]
        eb = end_bin[split]
        coef = coefficients[split]
        left = (edges[sb + 1] - starts[split]) * coef
        right = (ends[split] - edges[eb]) * coef
        accumulated += np.bincount(sb, weights=left, minlength=n_bins)
        accumulated += np.bincount(eb, weights=right, minlength=n_bins)

        has_middle = eb > sb + 1
        difference = np.zeros(n_bins + 1)
        np.add.at(difference, sb[has_middle] + 1, coef[has_middle])
        np.add.at(difference, eb[has_middle], -coef[has_middle])
        accumulated += np.cumsum(difference[:-1]) * np.diff(edges)

    return rate * accumulated


def expected_bin_counts(rate, edges, starts, ends, start_bin, end_bin, integrals):
    """Expected latent event allocations used by the likelihood gradient."""
    return back_project(rate, edges, starts, ends, start_bin, end_bin, 1.0 / integrals)


def project(values, rate, widths, edges, starts, ends, start_bin, end_bin):
    """Return the integral of rate * values over each observed interval."""
    weighted = rate * values
    cumulative = np.concatenate(([0.0], np.cumsum(weighted * widths)))
    start_value = cumulative[start_bin] + weighted[start_bin] * (starts - edges[start_bin])
    end_value = cumulative[end_bin] + weighted[end_bin] * (ends - edges[end_bin])
    return end_value - start_value


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


def newton_direction(theta, gradient, edges, starts, ends, smoothness, max_cg):
    """Approximately solve (-Hessian) d = gradient by preconditioned CG.

    Hessian-vector products are exact:
        H v = diag(allocated - n q) v - sum_i p_i (p_i . v) + n q (q . v) - smoothness L v
    where p_i is interval i's posterior allocation over bins and q the
    chromosome-wide share of each bin.  CG stops on negative curvature.
    """
    rate, widths, sb, eb, integrals = interval_integrals(theta, edges, starts, ends)
    n_events = len(starts)
    share = rate * widths / np.dot(rate, widths)
    allocated = back_project(rate, edges, starts, ends, sb, eb, 1.0 / integrals)
    likelihood_diagonal = allocated - n_events * share

    def negative_hessian(vector):
        projected = project(vector, rate, widths, edges, starts, ends, sb, eb) / integrals
        product = back_project(rate, edges, starts, ends, sb, eb, projected / integrals)
        product -= n_events * share * np.dot(share, vector)
        product -= likelihood_diagonal * vector
        steps = np.diff(vector)
        product[:-1] -= smoothness * steps
        product[1:] += smoothness * steps
        return product + vector.mean()  # pins the unidentifiable common log-rate

    # Jacobi preconditioner from an upper bound on sum_i p_ij^2.
    squared = back_project(rate, edges, starts, ends, sb, eb, 1.0 / integrals**2) * rate * widths
    inverse_diagonal = 1.0 / (np.maximum(squared, 1e-12) + 2.0 * smoothness)

    norm = np.linalg.norm(gradient)
    target = min(0.5, np.sqrt(norm)) * norm * 0.1
    solution = np.zeros_like(gradient)
    residual = gradient.copy()
    preconditioned = inverse_diagonal * residual
    search = preconditioned.copy()
    rz = np.dot(residual, preconditioned)
    for _ in range(max_cg):
        product = negative_hessian(search)
        curvature = np.dot(search, product)
        if curvature <= 0:
            break
        alpha = rz / curvature
        solution += alpha * search
        residual -= alpha * product
        if np.linalg.norm(residual) <= target:
            break
        preconditioned = inverse_diagonal * residual
        rz_next = np.dot(residual, preconditioned)
        search = preconditioned + (rz_next / rz) * search
        rz = rz_next
    return solution - solution.mean()


def fit_chromosome(starts, ends, chrom_length, bin_size, smoothness, iterations,
                   learning_rate, tolerance, rate_tolerance=1e-6, stability_window=3,
                   diagnostics=None):
    """Maximize the penalized likelihood with line-searched Newton-CG.

    Convergence requires max |gradient| <= tolerance (gradient of the log
    posterior with respect to each bin's log-rate, common-scale component
    removed) and a max relative change in normalized bin rates of at most
    rate_tolerance over the last stability_window iterations.  Relative
    objective change is never used: the likelihood carries a large additive
    baseline unrelated to fit stability.  Failure to converge raises.
    If diagnostics is a dict, it is filled with the final convergence summary.
    """
    if len(starts) == 0 or tolerance <= 0 or learning_rate <= 0 or rate_tolerance <= 0:
        raise ValueError("events, tolerances, and learning rate must be positive")
    if stability_window < 1:
        raise ValueError("stability window must be at least 1")
    edges = np.arange(0, chrom_length, bin_size, dtype=np.int64)
    edges = np.append(edges, chrom_length)
    widths = np.diff(edges)
    theta = np.zeros(len(edges) - 1)

    def log_share(values):
        # log of the normalized (scale-free) rate profile
        top = values.max()
        return values - top - np.log(np.dot(np.exp(values - top), widths))

    profile_history = [log_share(theta)]
    objective, gradient = objective_and_gradient(theta, edges, starts, ends, smoothness)
    gradient -= gradient.mean()  # remove the unidentifiable common log-rate
    summary = {}
    for iteration in range(iterations + 1):
        gradient_max = np.max(np.abs(gradient))
        window = profile_history[-(stability_window + 1):]
        rate_change = (np.max(np.abs(np.expm1(window[-1] - window[0])))
                       if len(window) == stability_window + 1 else np.inf)
        converged = gradient_max <= tolerance and rate_change <= rate_tolerance
        summary = dict(iterations=iteration, objective=objective,
                       max_gradient=gradient_max, max_rate_change=rate_change,
                       converged=converged)
        if converged or iteration == iterations:
            break
        direction = newton_direction(theta, gradient, edges, starts, ends,
                                     smoothness, max_cg=len(theta))
        slope = np.dot(gradient, direction)
        if slope <= 0 or not np.all(np.isfinite(direction)):
            direction = gradient / max(1.0, gradient_max)
            slope = np.dot(gradient, direction)
        step = learning_rate
        noise = 1e-10 * max(1.0, abs(objective))
        for _ in range(60):
            candidate = theta + step * direction
            candidate -= candidate.mean()
            try:
                with np.errstate(over="raise", invalid="raise", divide="raise"):
                    new_objective, new_gradient = objective_and_gradient(
                        candidate, edges, starts, ends, smoothness)
            except FloatingPointError:
                new_objective = -np.inf
            if np.isfinite(new_objective):
                new_gradient -= new_gradient.mean()
                if new_objective >= objective + 1e-4 * step * slope:
                    break
                # At the floating-point floor of the objective, accept a step
                # that still reduces the gradient.
                if (abs(new_objective - objective) <= noise and
                        np.max(np.abs(new_gradient)) < gradient_max):
                    break
            step *= 0.5
        else:
            raise RuntimeError(f"line search failed; max |gradient|={gradient_max:.3g}")
        theta, objective, gradient = candidate, new_objective, new_gradient
        profile_history.append(log_share(theta))
        profile_history = profile_history[-(stability_window + 1):]

    if diagnostics is not None:
        diagnostics.update(summary)
    if not summary["converged"]:
        raise RuntimeError(
            f"fit did not converge in {iterations} iterations: "
            f"max |gradient|={summary['max_gradient']:.3g} (tolerance {tolerance:g}), "
            f"max relative rate change over last {stability_window} iterations="
            f"{summary['max_rate_change']:.3g} (tolerance {rate_tolerance:g})")
    return edges, np.exp(theta), summary["iterations"], objective


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
    if args.tolerance <= 0 or args.rate_tolerance <= 0 or args.stability_window < 1:
        raise SystemExit("tolerances must be positive and the stability window at least 1")

    jri, targets, lengths = load_inputs(args)
    records = []
    summaries = []
    started = time.time()
    for chrom in sorted(jri["chr"].unique(), key=lambda value: int(value[3:])):
        if chrom not in targets or chrom not in lengths:
            raise SystemExit(f"ERROR: {chrom} has no Ogut target or v5 length; check chromosome names")
        subset = jri[jri["chr"] == chrom]
        starts = subset["start"].to_numpy(dtype=np.int64)
        ends = subset["end"].to_numpy(dtype=np.int64)
        chrom_length = int(lengths[chrom])
        if np.any(starts < 0) or np.any(ends > chrom_length):
            raise SystemExit(f"interval outside chromosome bounds on {chrom}")

        diagnostics = {}
        chrom_started = time.time()
        try:
            edges, relative_rate, iterations, objective = fit_chromosome(
                starts, ends, chrom_length, args.bin_size, args.smoothness,
                args.iterations, args.learning_rate, args.tolerance,
                args.rate_tolerance, args.stability_window, diagnostics,
            )
        except RuntimeError as error:
            print("\n" + "!" * 72, file=sys.stderr)
            print(f"ERROR: {chrom} NOT CONVERGED -- no map written.\n  {error}\n"
                  "  Raise --iterations or relax --tolerance/--rate-tolerance.",
                  file=sys.stderr)
            print("!" * 72, file=sys.stderr)
            raise SystemExit(1)
        summaries.append((chrom, diagnostics, time.time() - chrom_started))
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

    print(f"\nConvergence summary (tolerance {args.tolerance:g}, rate tolerance "
          f"{args.rate_tolerance:g} over {args.stability_window} iterations):",
          file=sys.stderr)
    print(f"  {'chrom':<6}{'iter':>6}{'max|grad|':>12}{'max rel dRate':>15}"
          f"{'converged':>11}{'sec':>7}", file=sys.stderr)
    for chrom, diag, seconds in summaries:
        print(f"  {chrom:<6}{diag['iterations']:>6}{diag['max_gradient']:>12.2e}"
              f"{diag['max_rate_change']:>15.2e}{str(diag['converged']):>11}"
              f"{seconds:>7.1f}", file=sys.stderr)
    print(f"  total fit time {time.time() - started:.1f} s", file=sys.stderr)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        for record in records:
            handle.write("\t".join(map(str, record)) + "\n")
    print(f"Wrote {len(records)} segments to {args.output}")


if __name__ == "__main__":
    main()
