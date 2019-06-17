import unittest
from bioseq import Blast
from bioseq import ProteinSeq
from bioseq.utils import *

db = "bioseq/resources/seqdump.txt"
fq = "bioseq/resources/source.fasta"
blosum = read_substitution_matrix_file("bioseq/resources/blosum62.mat")
sm = substitution_matrix("ATGC_", 1, -1)

class TestBlast(unittest.TestCase):
    def test_constructor(self):
        b = Blast(db, 3)
        self.assertIsInstance(b, Blast)
        self.assertEqual(b.k, 3)
        self.assertEqual(len(b.db), 21)

        b = Blast(list(read_fasta(db)), 3)
        self.assertIsInstance(b, Blast)
        self.assertEqual(b.k, 3)
        self.assertEqual(len(b.db), 21)

    def test_forward_matching(self):
        b = Blast(db, 3)
        q = next(read_fasta(fq))
        self.assertEqual(b.forward_matching(q, q, (0, 0), 0.5), (len(q) - 3, len(q) - 3))
        self.assertEqual(b.forward_matching(ProteinSeq("PFMI"), ProteinSeq("PFMI"), (0, 0), 0.5), (1, 1))
        self.assertEqual(b.forward_matching(ProteinSeq("PFMI"), ProteinSeq("PFMI"), (3, 3), 0.5), (0, 0))
        self.assertEqual(b.forward_matching(ProteinSeq("PFMIPPM"), ProteinSeq("PFMIPPF"), (0, 0), 0.5), (3, 4))
        self.assertEqual(b.forward_matching(ProteinSeq("PFMIPPM"), ProteinSeq("IIIIII"), (0, 0), 0.5), (1, 3))
        self.assertEqual(b.extend_hit(ProteinSeq("PFMI"), ProteinSeq("PFMI"), (0, 0), 0.5), 1.0)
        self.assertEqual(b.extend_hit(ProteinSeq("PFMI"), ProteinSeq("PFMI"), (3, 3), 0.5), 1.0)
        self.assertEqual(b.extend_hit(ProteinSeq("PFMIPPM"), ProteinSeq("PFMIPPF"), (0, 0), 0.5), 0.75)

    def test_search(self):
        b = Blast(db, 3)
        q = next(read_fasta(fq))
        similar = b.search(q, top=5)
        self.assertEqual(len(similar), 5)
        self.assertEqual(similar[0][0].name, "XP_016791527.1 cytochrome P450 2B6 [Pan troglodytes]")
        self.assertAlmostEqual(similar[0][1], 0.9901960784313726, places=5)
        self.assertAlmostEqual(similar[1][1], 0.9901960784313726, places=5)
        self.assertAlmostEqual(similar[2][1], 0.9852941176470589, places=5)

    #     m = MSA(blosum, -8)
    #     seqs = [ProteinSeq("PHWASW"), ProteinSeq("HPHWA")]
    #     self.assertEqual(str(m.align(seqs)), "(HPHWASW, [_PHWASW, HPHWA__])")



if __name__ == '__main__':
    unittest.main()
