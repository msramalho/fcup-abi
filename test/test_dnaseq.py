import unittest
from bioseq import DNASeq, RNASeq, ProteinSeq


class TestDNASeq(unittest.TestCase):
    def test_constructor(self):
        s = DNASeq("ATATATCGCG")
        self.assertIsInstance(s, DNASeq)
        self.assertEqual("ATATATCGCG", s.sequence)
        self.assertEqual("DNA", s.seq_type)

    def test_custom_assert_valid(self):
        s = DNASeq("ATATATCGCG")
        self.assertTrue(s._assert_valid_sequence())
        s.sequence += "X"
        self.assertFalse(s._assert_valid_sequence())

    def test_reverse_complement(self):
        s = DNASeq("ATCG")
        self.assertEqual("CGAT", s.reverse_complement())

    def test_transcription(self):
        dna = DNASeq("ATCG")
        rna = dna.transcription()
        self.assertEqual("AUCG", rna.sequence)
        self.assertEqual("RNA", rna.seq_type)
        self.assertIsInstance(rna, RNASeq)

    def test_translation(self):
        dna = DNASeq("GCTGCCTGTATTTAG")
        pro = dna.translation()
        self.assertEqual("AACI_", pro.sequence)
        self.assertEqual("PROTEIN", pro.seq_type)
        self.assertIsInstance(pro, ProteinSeq)
        self.assertEqual("AACI_", dna.translation(0).sequence)  # expecting the same as the default
        self.assertEqual("ACI_", dna.translation(3).sequence)  # skip the first codon
        dna = DNASeq("AGCTGCCTGTATTTAG")
        self.assertEqual("AACI_", dna.translation(1).sequence) # test different start
        dna = DNASeq("AAGCTGCCTGTATTTAG")
        self.assertEqual("AACI_", dna.translation(2).sequence) # test different start

    def test_codon_usage(self):
        dna = DNASeq("GCTGCTGCTGCCGCCGCAGCAGCGGCGGCG")
        self.assertDictEqual({'GCA': 0.2, 'GCC': 0.2, 'GCG': 0.3, 'GCU': 0.3}, dna.codon_usage("A"))

if __name__ == '__main__':
    unittest.main()
