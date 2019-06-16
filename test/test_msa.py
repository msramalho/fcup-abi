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
        m = MSA(sm2, -1)

        seqs = [DNASeq("ATAGC"), DNASeq("AACC")]
        c, aligned = m.align(seqs)
        self.assertEqual("ATACC", str(c))
        self.assertEqual(len(aligned), 2)

        seqs = [DNASeq("ATAGC"), DNASeq("AACC"), DNASeq("ATGAC")]
        c, aligned = m.align(seqs)
        self.assertEqual("ATGACC", str(c))
        self.assertEqual(len(aligned), 3)

    def test_consensus(self):
        m = MSA(sm2, -1)
        c = m.consensus([DNASeq("ATAGC"), DNASeq("A_ACC")], DNASeq)
        self.assertEqual(type(c), bioseq.dnaseq.DNASeq)
        self.assertEqual(str(c), "ATACC")

    def test_get_col(self):
        m = MSA(sm2, -1)
        seqs = [DNASeq("ATAGC"), DNASeq("A_ACC")]
        self.assertEqual(m.get_col(seqs, 0), ["A", "A"])
        self.assertEqual(m.get_col(seqs, 1), ["T", "_"])
        

if __name__ == '__main__':
    unittest.main()
