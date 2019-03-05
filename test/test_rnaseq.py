import unittest
from bioseq import RNASeq, ProteinSeq


class TestRNASeq(unittest.TestCase):
    def test_constructor(self):
        s = RNASeq("AUGCGAU")
        self.assertIsInstance(s, RNASeq)
        self.assertEqual("AUGCGAU", s.sequence)
        self.assertEqual("RNA", s.seq_type)

    def test_custom_assert_valid(self):
        s = RNASeq("AUGCGAU")
        self.assertTrue(s._assert_valid_sequence())
        s.sequence += "T"
        self.assertFalse(s._assert_valid_sequence())

    def test_reverse_complement(self):
        s = RNASeq("AUCG")
        self.assertEqual("CGAU", s.reverse_complement())

    def test_translation(self):
        rna = RNASeq("GCUGCCUGUAUUUAG")
        pro = rna.translation()
        self.assertEqual("AACI_", pro.sequence)
        self.assertEqual("PROTEIN", pro.seq_type)
        self.assertIsInstance(pro, ProteinSeq)
        self.assertEqual("AACI_", rna.translation(0).sequence)  # expecting the same as the default
        self.assertEqual("ACI_", rna.translation(3).sequence)  # skip the first codon
        rna = RNASeq("CGCUGCCUGUAUUUAG")
        self.assertEqual("AACI_", rna.translation(1).sequence)  # test different start
        rna = RNASeq("CCGCUGCCUGUAUUUAG")
        self.assertEqual("AACI_", rna.translation(2).sequence)  # test different start

    def test_codon_usage(self):
        rna = RNASeq("GCUGCUGCUGCCGCCGCAGCAGCGGCGGCG")
        self.assertDictEqual({'GCA': 0.2, 'GCC': 0.2, 'GCG': 0.3, 'GCU': 0.3}, rna.codon_usage("A"))


if __name__ == '__main__':
    unittest.main()
