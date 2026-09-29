"""Regression checks for haploid VCF parsing and BED membership."""

import io
import unittest

import numpy as np

from haploid_pi import accumulate_stream


class HaploidPiTests(unittest.TestCase):
    def accumulate(self, records, chunk_bytes=100000, restrict=None):
        header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n"
        stream = io.BytesIO((header + records).encode())
        return accumulate_stream(stream, "fixture", {"p": ["a", "b"]},
                                 {"chr1": (np.array([0, 3]), np.array([3, 10]))},
                                 chunk_bytes, False, restrict)[1]["p"]

    def record(self, pos=1, fmt="GT:DP", a="0:8", b="0:9", chrom="chr1", alt="C"):
        return f"{chrom}\t{pos}\t.\tA\t{alt}\t.\t.\t.\t{fmt}\t{a}\t{b}\n"

    def test_format_fields_do_not_become_alleles(self):
        acc = self.accumulate(self.record())
        np.testing.assert_array_equal(acc["diffs"], [0, 0])
        np.testing.assert_array_equal(acc["pairs"], [1, 0])

    def test_gt_order_missing_and_variable_width_windows(self):
        records = (self.record(3, "DP:GT", "8:0", "9:1")
                   + self.record(4, "GT", "0", ".")
                   + self.record(10, "GT", "0", "1").rstrip("\n"))
        acc = self.accumulate(records)
        np.testing.assert_array_equal(acc["diffs"], [1, 1])
        np.testing.assert_array_equal(acc["pairs"], [1, 1])

    def test_multidigit_alleles(self):
        acc = self.accumulate(self.record(fmt="GT", a="10", b="0", alt=",".join(["C"] * 10)))
        self.assertEqual(acc["diffs"][0], 1)

    def test_reject_invalid_records_even_outside_restriction(self):
        invalid = [self.record(pos=0), self.record(pos=11),
                   self.record(a="0/1:8"), self.record(a="2:8"),
                   self.record(fmt="DP", a="8", b="9")]
        for record in invalid:
            with self.subTest(record=record), self.assertRaises(SystemExit):
                self.accumulate(record, restrict={"chr1": np.array([2])})

    def test_checks_chromosome_within_and_between_chunks(self):
        for size in [1, 100000]:
            with self.subTest(size=size), self.assertRaises(SystemExit):
                self.accumulate(self.record() + self.record(chrom="chr2"), size)


if __name__ == "__main__":
    unittest.main()
