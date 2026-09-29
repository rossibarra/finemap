"""Regression checks for observation-supported crossover intervals and QC."""
import unittest
from unittest.mock import patch

import hmm_co_pipeline as hmm


class HMMTests(unittest.TestCase):
    def markers(self, observations):
        return [hmm.Marker(1, i * 1_000_000, str(i), gt)
                for i, gt in enumerate(observations)]

    def test_missing_calls_do_not_narrow_interval_and_ids_survive(self):
        events = hmm.build_hmm_events("map", self.markers("AAA---BBB"), [6])
        self.assertEqual(len(events), 1)
        self.assertEqual((events[0].left_coordinate, events[0].right_coordinate),
                         (2_000_000, 6_000_000))
        self.assertEqual(events[0].sample_id, "map_ind007")

    def test_discordant_calls_and_unsupported_runs(self):
        self.assertEqual(list(hmm.supported_breakpoints(
            self.markers("AABABB"), [0, 0, 0, 1, 1, 1], 0)), [(1, 4)])
        self.assertEqual(list(hmm.supported_breakpoints(
            self.markers("AA--AA"), [0, 0, 1, 1, 0, 0], 0)), [])

    def test_chromosome_boundary_is_not_a_neighbor(self):
        markers = self.markers("ABA")
        markers[-1].chromosome = 2
        self.assertEqual(hmm.isolated_flip_rate(markers, 1), 0)
        markers[-1].chromosome = 1
        self.assertEqual(hmm.isolated_flip_rate(markers, 1), 1)

    def test_flip_qc_uses_retained_samples_only(self):
        markers = self.markers(["ABA", "ABB", "ABA", "AB-", "AB-", "AB-"])
        with patch.object(hmm, "isolated_flip_rate", wraps=hmm.isolated_flip_rate) as rate:
            cleaned, kept, summary = hmm.clean_population(markers)
        self.assertEqual(kept, [0, 1])
        self.assertEqual(summary["removed_markers_isolated_flip"], 0)
        self.assertTrue(all(len(m.genotypes) == 2 for m in rate.call_args.args[0]))
        self.assertEqual(len(cleaned), 6)

    def test_ragged_genotypes_count_as_missing(self):
        self.assertEqual(hmm.marker_missingness("A", [0, 1]), 0.5)
        self.assertEqual(hmm.allele_balance("A", [0, 1]), 1)


if __name__ == "__main__":
    unittest.main()
