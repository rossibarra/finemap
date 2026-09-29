"""Regression tests for fit stationarity and bootstrap genomic gaps."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from build_hierarchical_finemap import fit_chromosome, objective_and_gradient
from pi_vs_recombination import physical_blocks, block_bootstrap_ci


class HierarchicalFitTests(unittest.TestCase):
    def setUp(self):
        self.starts = np.array([0, 10, 20, 25, 70, 80])
        self.ends = np.array([10, 20, 50, 100, 95, 100])

    def test_gradient_matches_finite_difference(self):
        edges = np.arange(0, 101, 10)
        theta = np.linspace(-0.4, 0.7, 10)
        _, gradient = objective_and_gradient(theta, edges, self.starts, self.ends, 2)
        differences = []
        for i in range(10):
            delta = np.eye(10)[i] * 1e-5
            lo, _ = objective_and_gradient(theta - delta, edges, self.starts, self.ends, 2)
            hi, _ = objective_and_gradient(theta + delta, edges, self.starts, self.ends, 2)
            differences.append((hi - lo) / 2e-5)
        np.testing.assert_allclose(gradient, differences, atol=1e-8)

    def test_fit_is_stationary_and_stable(self):
        fits = []
        for tolerance in (1e-5, 1e-7):
            edges, rate, _, _ = fit_chromosome(
                self.starts, self.ends, 100, 10, 2, 1000, 1, tolerance)
            _, gradient = objective_and_gradient(np.log(rate), edges,
                                                  self.starts, self.ends, 2)
            self.assertLessEqual(np.max(np.abs(gradient)), tolerance)
            fits.append(rate / np.dot(rate, np.diff(edges)))
        np.testing.assert_allclose(*fits, rtol=1e-4)

    def test_unfinished_fit_fails(self):
        with self.assertRaisesRegex(RuntimeError, "did not converge"):
            fit_chromosome(self.starts, self.ends, 100, 10, 2, 1, 1, 1e-7)


class BootstrapTests(unittest.TestCase):
    def test_gaps_and_chromosomes_do_not_join(self):
        df = pd.DataFrame({"chrom": ["chr1"] * 4 + ["chr2"],
                           "start": [0, 100, 5100, 9900, 0],
                           "end": [100, 200, 5200, 10000, 100],
                           "x": [1, 2, 3, 4, 5], "y": [2, 1, 4, 3, 6]})
        blocks = physical_blocks(df, 5000)
        self.assertEqual([b.start.tolist() for b in blocks], [[0, 100], [5100, 9900], [0]])
        # Imported clients can still supply window counts, with identical results.
        legacy = block_bootstrap_ci(df, "x", "y", 50, 100, 1)
        explicit = block_bootstrap_ci(df, "x", "y", None, 100, 1, block_bp=5000)
        self.assertEqual(legacy, explicit)


if __name__ == "__main__":
    unittest.main()
