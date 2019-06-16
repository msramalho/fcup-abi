import unittest
from bioseq import Blast, ProteinSeq, DNASeq, MSA
from bioseq.utils import *

db = "bioseq/resources/seqdump.txt"
fq = "bioseq/resources/source.fasta"
sm = read_substitution_matrix_file("test/blosum62.mat")
sm2 = substitution_matrix("ATGC", 1, -1)


class TestMSA(unittest.TestCase):
    def test_constructor(self):
        m = MSA(sm, -3)
        self.assertEqual(m.sm, sm)
        self.assertEqual(m.g, -3)

    def test_align(self):
        seqs = [DNASeq("ATAGC"), DNASeq("AACC")]
        m = MSA(sm2, -1)
        self.assertEqual("ATACC", str(m.align(seqs)))
        seqs = [DNASeq("ATAGC"), DNASeq("AACC"), DNASeq("ATGAC")]
        self.assertEqual("ATGACC", str(m.align(seqs)))

    def test_consensus(self):
        m = MSA(sm2, -1)
        c = m.consensus(DNASeq("ATAGC"), DNASeq("A_ACC"))
        self.assertEqual(type(c), bioseq.dnaseq.DNASeq)
        self.assertEqual(str(c), "ATACC")
        
    def test_consensus(self):
        m = MSA(sm2, -1)
        self.assertEqual(m.consensus_col("A", GAP), "A")
        self.assertEqual(m.consensus_col(GAP, "A"), "A")
        self.assertEqual(m.consensus_col("A", "A"), "A")
        self.assertEqual(m.consensus_col("A", "C"), "A")
        self.assertEqual(m.consensus_col("C", "A"), "A")
        self.assertEqual(m.consensus_col("C", "C"), "C")
        self.assertEqual(m.consensus_col("C", "G"), "C")
        self.assertEqual(m.consensus_col("G", "C"), "C")
        self.assertEqual(m.consensus_col("G", "G"), "G")
        self.assertEqual(m.consensus_col("G", "T"), "G")
        self.assertEqual(m.consensus_col("T", "G"), "G")
        self.assertEqual(m.consensus_col("T", "T"), "T")
        self.assertEqual(m.consensus_col("C", GAP), "C")
        self.assertEqual(m.consensus_col(GAP, "C"), "C")
        self.assertEqual(m.consensus_col("G", GAP), "G")
        self.assertEqual(m.consensus_col(GAP, "G"), "G")
        self.assertEqual(m.consensus_col("T", GAP), "T")
        self.assertEqual(m.consensus_col(GAP, "T"), "T")


if __name__ == '__main__':
    unittest.main()
