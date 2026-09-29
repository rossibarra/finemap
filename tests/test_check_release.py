"""Synthetic-fixture tests for scripts/check_release.py."""
import io
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import check_release as cr

LENGTH = 10_000
# Per chromosome: one interval per source; LRv4 and TEOv4 overlap exactly.
INTERVALS = [('RMv2', 100, 2000), ('EUROv2', 2000, 5000), ('LRv4', 6000, 9000),
             ('TEOv4', 6000, 9000)]


def ogut_total(i):
    return 10.0 + i


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def make_fixture(root):
    """Write a small, fully consistent release under root."""
    root = Path(root)
    write(root / cr.FAI, ''.join(f'chr{i}\t{LENGTH}\t0\t80\t81\n' for i in range(1, 11))
          + 'scaf_1\t500\t0\t80\t81\n')
    write(root / cr.OGUT, 'SNP_ID,SNP_newID,chromosome,position,cM\n' + ''.join(
        f'S{i}a,M{i}a,{i},10,-1.0\nS{i}b,M{i}b,{i},900,{ogut_total(i) - 1.0}\n'
        for i in range(1, 11)))

    # HMM: one event per chromosome for sample S{i}; EUROv2_{i} is row i.
    write(root / cr.HMM_EVENTS,
          'sample_id\tmap\tchromosome\tleft_coordinate\tright_coordinate\n' + ''.join(
              f'S{i}\tM\t{i}\t2000\t5000\n' for i in range(1, 11)))
    write(root / cr.HMM_ROSTER, 'sample_id\tmap\tchromosome\n' + ''.join(
        f'S{i}\tM\t{i}\n' for i in range(1, 11)))

    jri, audit, per_source = [], [], {s: [] for s in cr.SOURCES}
    audit.append('source\tid\tsample\tsource_chrom\tsource_start_1based\tsource_end_1based\t'
                 'left\tright\tstatus\tlength_ratio\textreme_length_change')
    for i in range(1, 11):
        for src, s, e in INTERVALS:
            key = f'{src}_{i:06d}'
            sample = f'S{i}' if src == 'EUROv2' else f'T{i}'
            jri.append(f'Chr{i}\t{s}\t{e}\t{sample}\t{key}')
            per_source[src].append(f'Chr{i}\t{s}\t{e}\t{sample}')
            tag = 'v2' if src in ('RMv2', 'EUROv2') else src
            audit_id = key if tag == 'v2' else f'{i}_{s}_{e}_{sample}'
            audit.append(f'{tag}\t{audit_id}\t{sample}\t{i}\t{s}\t{e}\t\t\tkept\t1.0\t0')
    audit.append('v2\tRMv2_999999\tX\t1\t5\t50\t\t\tunmapped_endpoint\t\t0')
    write(root / cr.JRI, '\n'.join(jri) + '\n')
    write(root / cr.AUDIT, '\n'.join(audit) + '\n')
    write(root / cr.JRI_RM_EURO, ''.join(
        f'{r}\t{s}_{n}\n' for s in ('RMv2', 'EUROv2') for n, r in enumerate(per_source[s])))
    for src, rel in cr.XO_V5.items():
        write(root / rel, '\n'.join(per_source[src]) + '\n')

    # finemap_v5: sweep of the intervals above, scaled to the Ogut total.
    fm, hier = [], []
    for i in range(1, 11):
        segs = [(100, 2000, 1 / 1900), (2000, 5000, 1 / 3000), (6000, 9000, 2 / 3000)]
        scale = ogut_total(i) / sum((e - s) * w for s, e, w in segs)
        cm = 0.0
        hap = [(0, 0.0)]
        for s, e, w in segs:
            end_cm = cm + (e - s) * w * scale
            fm.append(f'Chr{i}\t{s}\t{e}\t{cm!r}\t{end_cm!r}\t{w * scale * 1e6!r}')
            hap += [(s, cm), (e, end_cm)]
            cm = end_cm
        hap.append((LENGTH, cm))
        points = dict(hap)
        pos = sorted(points)
        lines = ['Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)']
        for j, p in enumerate(pos):
            rate = 0.0 if j == len(pos) - 1 else \
                (points[pos[j + 1]] - points[p]) / (pos[j + 1] - p) * 1e6
            lines.append(f'chr{i}\t{p}\t{rate:.10g}\t{points[p]:.10g}')
        write(root / f'data/hapmap/chr{i}.hapmap.tsv', '\n'.join(lines) + '\n')
        rate = ogut_total(i) / LENGTH
        for b in range(0, LENGTH, 2500):
            hier.append(f'Chr{i}\t{b}\t{b + 2500}\t{rate * b!r}\t{rate * (b + 2500)!r}\t'
                        f'{rate * 1e6!r}')
    write(root / cr.FINEMAP, '\n'.join(fm) + '\n')
    write(root / cr.HIER, '\n'.join(hier) + '\n')
    write(root / cr.HOT30, ''.join(f'chr{i}\t2100\t2200\t0.1\t40\t2\n' for i in range(1, 11)))
    write(root / cr.HOT1KB, ''.join(f'chr{i}\t6000\t7000\t25\n' for i in range(1, 11)))


class ReleaseCheckTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='finemap-release-')
        self.root = Path(self.tmp.name)
        make_fixture(self.root)

    def tearDown(self):
        self.tmp.cleanup()

    def run_check(self, *extra):
        out = io.StringIO()
        code = cr.main(['--root', str(self.root), *extra], out=out)
        return code, out.getvalue()

    def edit(self, rel, old, new, count=1):
        path = self.root / rel
        text = path.read_text()
        self.assertIn(old, text)
        path.write_text(text.replace(old, new, count))

    def status(self, output, name):
        for line in output.splitlines():
            if line[6:].startswith(name):
                return line[:4]
        self.fail(f'no line for {name}:\n{output}')

    def test_consistent_fixture_passes(self):
        code, out = self.run_check()
        self.assertEqual(code, 0, out)
        self.assertIn('WARN  provenance manifest', out)
        self.assertIn('RELEASE CHECK PASSED', out)

    def test_decreasing_cm_fails(self):
        lines = (self.root / cr.FINEMAP).read_text().splitlines()
        f = lines[1].split('\t')
        f[4] = repr(float(f[3]) - 0.5)  # cM_end below cM_start
        lines[1] = '\t'.join(f)
        (self.root / cr.FINEMAP).write_text('\n'.join(lines) + '\n')
        code, out = self.run_check()
        self.assertEqual(code, 1)
        self.assertEqual(self.status(out, 'finemap_v5.bed'), 'FAIL')
        self.assertIn('decreasing cM', out)

    def test_ogut_total_mismatch_fails_both_maps(self):
        self.edit(cr.OGUT, 'S3b,M3b,3,900,12.0', 'S3b,M3b,3,900,12.5')
        code, out = self.run_check()
        self.assertEqual(code, 1)
        self.assertEqual(self.status(out, 'finemap_v5.bed'), 'FAIL')
        self.assertEqual(self.status(out, 'finemap_hierarchical_v5.bed'), 'FAIL')
        self.assertIn('Chr3: final cM', out)
        self.assertEqual(self.status(out, 'hapmap exports'), 'PASS')

    def test_audit_kept_count_mismatch_fails(self):
        self.edit(cr.AUDIT, 'LRv4\t2_6000_9000_T2\tT2\t2\t6000\t9000\t\t\tkept',
                  'LRv4\t2_6000_9000_T2\tT2\t2\t6000\t9000\t\t\tdifferent_strand')
        code, out = self.run_check()
        self.assertEqual(code, 1)
        self.assertEqual(self.status(out, 'lift-over audit'), 'FAIL')
        self.assertIn('LRv4: 9 kept in audit vs 10 LRv4_ IDs', out)

    def test_hmm_rerun_detected(self):
        # Extra HMM event at the top: EUROv2 raw count and row identities no longer match.
        self.edit(cr.HMM_EVENTS, 'right_coordinate\n', 'right_coordinate\nS1\tM\t1\t10\t20\n')
        code, out = self.run_check()
        self.assertEqual(code, 1)
        self.assertEqual(self.status(out, 'jri_v5.bed'), 'FAIL')
        self.assertIn('EUROv2: 10 audit rows vs 11 events', out)

    def test_hapmap_disagreement_fails(self):
        self.edit('data/hapmap/chr5.hapmap.tsv', 'chr5\t2000\t', 'chr5\t2001\t')
        code, out = self.run_check()
        self.assertEqual(self.status(out, 'hapmap exports'), 'FAIL')
        self.assertIn('map boundaries absent', out)

    def test_hotspot_outside_map_fails(self):
        self.edit(cr.HOT1KB, 'chr4\t6000\t7000', 'chr4\t8500\t9500')
        code, out = self.run_check()
        self.assertEqual(self.status(out, 'hotspots_1kb_20x_v5.bed'), 'FAIL')

    def test_stale_manifest_detected(self):
        code, out = self.run_check('--write-manifest')
        self.assertEqual(code, 0, out)
        code, out = self.run_check()
        self.assertEqual(code, 0, out)
        self.assertIn('no stale products', out)

        # Rewrite jri_v5.bed (consistent content change) without rebuilding the maps.
        self.edit(cr.JRI, 'Chr1\t100\t2000\tT1', 'Chr1\t100\t2000\tT1b')
        code, out = self.run_check()
        self.assertEqual(code, 1)
        self.assertEqual(self.status(out, 'dependency order'), 'FAIL')
        for rel in (cr.FINEMAP, cr.HIER, cr.HOT30, cr.HOT1KB, 'data/hapmap/chr1.hapmap.tsv'):
            self.assertIn(f'{rel} is stale', out)
        self.assertNotIn(f'{cr.HMM_EVENTS} is stale', out)
        self.assertNotIn(f'{cr.JRI} is stale', out)

    def test_manifest_refused_when_inconsistent(self):
        (self.root / cr.AUDIT).unlink()
        code, out = self.run_check('--write-manifest')
        self.assertEqual(code, 1)
        self.assertFalse((self.root / cr.MANIFEST).exists())
        self.assertIn('results/liftover_audit.tsv missing', out)


if __name__ == '__main__':
    unittest.main()
