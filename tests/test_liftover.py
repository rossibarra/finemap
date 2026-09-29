"""Coordinate and structural-boundary regressions for crossover lifting."""
import csv
import io
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from build_jri_v5 import lift_endpoints


class LiftTests(unittest.TestCase):
    def run_lift(self, chain, rows):
        with tempfile.TemporaryDirectory(prefix='finemap-lift-') as directory:
            work = Path(directory)
            path = work / 'test.chain'
            path.write_text(chain)
            audit = io.StringIO()
            kept = lift_endpoints(rows, path, work, 'test', csv.writer(audit, delimiter='\t'))
            return kept, audit.getvalue()

    def test_forward_and_reverse_one_based_markers(self):
        chain = ('chain 1 1 100 + 0 20 chr1 1000 + 100 120 1\n20\n\n'
                 'chain 1 1 100 + 20 40 chr1 1000 - 200 220 2\n20\n')
        kept, audit = self.run_lift(chain, [('1', 1, 20, 's', 'f'), ('1', 21, 40, 's', 'r')])
        self.assertEqual(kept, {'f': (1, 100, 120), 'r': (1, 780, 800)})
        self.assertIn('chr1:799:-:2', audit)

    def test_different_chains_and_orientations_rejected(self):
        chain = ('chain 1 1 100 + 0 20 chr1 1000 + 100 120 1\n20\n\n'
                 'chain 1 1 100 + 20 40 chr1 1000 + 200 220 2\n20\n\n'
                 'chain 1 1 100 + 40 60 chr1 1000 - 300 320 3\n20\n')
        kept, audit = self.run_lift(chain, [('1', 1, 40, 's', 'chain'), ('1', 1, 60, 's', 'strand')])
        self.assertEqual(kept, {})
        self.assertIn('different_chain', audit)
        self.assertIn('different_strand', audit)

    def test_two_hits_on_one_endpoint_do_not_replace_missing_other(self):
        chain = ('chain 1 1 100 + 0 20 chr1 1000 + 100 120 1\n20\n\n'
                 'chain 1 1 100 + 0 20 chr1 1000 + 200 220 2\n20\n')
        kept, audit = self.run_lift(chain, [('1', 1, 40, 's', 'ambiguous')])
        self.assertFalse(kept)
        self.assertIn('ambiguous_endpoint', audit)

    def test_same_chain_gap_retained_and_expansion_flagged(self):
        chain = 'chain 1 1 100 + 0 21 chr1 1000 + 0 520 1\n10 1 500\n10\n'
        kept, audit = self.run_lift(chain, [('1', 1, 21, 's', 'gap')])
        self.assertEqual(kept['gap'], (1, 0, 520))
        self.assertTrue(audit.rstrip().endswith('\t1'))

    def test_block_end_is_exclusive(self):
        chain = 'chain 1 1 100 + 0 20 chr1 1000 + 100 120 1\n20\n'
        kept, audit = self.run_lift(chain, [('1', 20, 21, 's', 'gap')])
        self.assertFalse(kept)
        self.assertIn('unmapped_endpoint', audit)


if __name__ == '__main__':
    unittest.main()
