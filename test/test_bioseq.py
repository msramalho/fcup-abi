import io
import os
import unittest
import unittest.mock
from bioseq import BioSeq
from bioseq.utils import *


class TestBioSeq(unittest.TestCase):
    def test_constructor(self):
        s = BioSeq("ATATat", "DNa")
        self.assertIsInstance(s, BioSeq)
        self.assertEqual("ATATAT", s.sequence)
        self.assertEqual("DNA", s.seq_type)
        self.assertRaises(Exception, BioSeq, "ATAT", "smth")

    def test_frequency(self):
        s = BioSeq("AACCCCCTGG", "DNA")
        self.assertDictEqual({"A": 0.2, "C": 0.5, "T": 0.1, "G": 0.2}, s.frequency())

    def test_gc_content(self):
        s = BioSeq("AACCCCCTGG", "DNA")
        self.assertEqual(0.70, s.gc_content())
        s.sequence = "AA"
        self.assertEqual(0.0, s.gc_content())
        s.sequence = "GC"
        self.assertEqual(1.0, s.gc_content())

    def test_assert_seq_type(self):
        s = BioSeq("ATATAT", "DNA")
        s = BioSeq("ATATAT", "RNA")
        s = BioSeq("ATATAT", "PROTEIN")
        self.assertRaises(Exception, BioSeq, "ATATA", "asd")

    def test_add_gap(self):
        s = BioSeq("ATGC", "DNA")
        self.assertEqual("ATGC", str(s))
        s.add_gap(2)
        self.assertEqual("AT_GC", str(s))


    def test_eq(self):
        s1 = BioSeq("ATATAT", "DNA")
        s2 = BioSeq("ATATAT", "DNA")
        self.assertTrue(s1==s2)
        self.assertFalse(s1!=s2)
        s3 = BioSeq("ATATATX", "DNA")
        s4 = BioSeq("ATATAT", "PROTEIN")
        self.assertFalse(s1==s3)

    def test_str(self):
        s = BioSeq("ATATAT", "DNA")
        self.assertEqual("ATATAT", str(s))

    def test_repr(self):
        s = [BioSeq("ATATAT", "DNA")]
        self.assertEqual("[ATATAT]", str(s))

    def test_getitem(self):
        s = BioSeq("ATATAT", "DNA")
        self.assertEqual("A", s[0])
        self.assertEqual("T", s[-1])

    def test_getslice(self):
        s = BioSeq("ATATAT", "DNA")
        self.assertEqual("A", s[0:1])
        self.assertEqual("T", s[-1:])
        self.assertEqual("AT", s[0:2])
        self.assertEqual("A", s[0:1:1])
        self.assertEqual("T", s[-1::1])
        self.assertEqual("AAA", s[0::2])
        # explicit invocation of function
        self.assertEqual("A", s.__getslice__(0, 1))
        self.assertEqual("T", s.__getslice__(-1, len(s)))
        self.assertEqual("AT", s.__getslice__(0, 2))
        self.assertEqual("A", s.__getslice__(0, 1, 1))
        self.assertEqual("T", s.__getslice__(-1, len(s), 1))
        self.assertEqual("AAA", s.__getslice__(0, len(s), 2))

    def test_len(self):
        s = BioSeq("ATATAT", "DNA")
        self.assertEqual(6, len(s))

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def _test_pretty_print(self, bioseq, mock_stdout):
        bioseq.pretty_print()
        output = mock_stdout.getvalue()
        self.assertEqual(len(output), 186)
        self.assertIn(bioseq.sequence, output)

    def test_pretty_print(self):
        s = BioSeq("ATATACAGATGAT", "DNA")
        self._test_pretty_print(s)

    def test_save_load(self):
        filename = "test_save_load.tmp"
        s = BioSeq("ATATACAGATGAT", "DNA")
        s.save(filename)
        s2 = BioSeq("", "RNA")
        s2.load(filename)
        self.assertEqual(s.seq_type, s2.seq_type)
        self.assertEqual(s.sequence, s2.sequence)
        os.remove(filename)

    def test_hamming_distance(self):
        s1 = BioSeq("ATATACAGATGAT", "DNA")
        s2 = BioSeq("ATATACAGATGAT", "DNA")
        self.assertEqual(0, s1.hamming_distance(s2))
        self.assertEqual(0, s2.hamming_distance(s1))
        s3 = BioSeq("ATATACAGATGAX", "DNA")
        self.assertEqual(1, s1.hamming_distance(s3))
        self.assertEqual(1, s3.hamming_distance(s1))
        s4 = BioSeq("AAA", "DNA")
        self.assertRaises(Exception, s1.hamming_distance, s4)

    def test_dot_plot(self):
        s1 = BioSeq("ATAT", "DNA")
        s2 = BioSeq("ATAT", "DNA")
        m = s1.dot_plot(s2)
        self.assertEqual(m.sum(), 8)
        s3 = BioSeq("XXXT", "DNA")
        m = s1.dot_plot(s3)
        self.assertEqual(m.sum(), 2)

    def test_score_seq(self):
        sm = substitution_matrix("ATCG", 2, -3)
        s1 = BioSeq("ATAT", "DNA")
        s2 = BioSeq("ATAG", "DNA")
        self.assertEqual(3, s1.score_seq(s2, sm, 3))

    def test_score_affine_gap(self):
        sm = read_substitution_matrix_file("test/blosum62.mat")
        s1 = BioSeq("LGPSSGCASRIWTKSA", "PROTEIN")
        s2 = BioSeq("TGPS_G__S_IWSKSG", "PROTEIN")
        self.assertEqual(33, s1.score_affine_gap(s2, sm, -8, -2))

    def test_global_align(self):
        # example from page 39 from slides
        s1 = BioSeq("PHSWG", "PROTEIN")
        s2 = BioSeq("HGWAG", "PROTEIN")
        sm = read_substitution_matrix_file("test/blosum62.mat")
        s, t = s1.global_align_multiple_solutions(s2, sm, -8)
        self.assertListEqual([-40, -24, -10, 3, 11, 9], s[-1])
        recover = list(s1.recover_global_align_multiple_solutions(s2, t))
        self.assertListEqual([('PHSW_G', '_HGWAG')], recover)


    def test_global_align2(self):
        # example from page 46 from slides
        s1 = BioSeq("PHSWG", "PROTEIN")
        s2 = BioSeq("HGWAG", "PROTEIN")
        sm = read_substitution_matrix_file("test/blosum62.mat")
        s, t = s1.local_align_multiple_solutions(s2, sm, -8)
        self.assertListEqual([0, 0, 6, 11, 19, 17], s[-1])
        recover = list(s1.recover_local_align_multiple_solutions(s2, t, s))
        self.assertListEqual([('HSW', 'HGW'), ('HSWG', 'HGWA')], recover)
        self.assertEqual(len(BioSeq.compare_pairwise_global_align([s1,s2], sm, -8)), 2)
        self.assertEqual(len(BioSeq.compare_pairwise_local_align([s1,s2], sm, -8)), 2)



if __name__ == '__main__':
    unittest.main(verbosity=2)
