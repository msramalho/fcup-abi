import unittest
from bioseq import *

blosum = read_substitution_matrix_file("bioseq/resources/blosum62.mat")
sm = substitution_matrix("ATGC_", 1, -1)


class TestPTree(unittest.TestCase):
    def test_constructor(self):
        p = PTree(blosum, -8)
        self.assertEqual(p.sm, blosum)
        self.assertEqual(p.g, -8)
        self.assertEqual(p.ignore, 0)

    def test_distance_matrix(self):
        p = PTree(sm, -1)
        seqs = [DNASeq("A_CATATC_AT_"), DNASeq("A_GATATT_AG_"), DNASeq("AACAGATC_T__"), DNASeq("G_CAT__CGATT")]
        m = p.distance_matrix(seqs)
        self.assertEqual(m[0], [0, 3, 3, 5])
        self.assertEqual(m[1], [0, 0, 6, 8])
        self.assertEqual(m[2], [0, 0, 0, 5])
        self.assertEqual(m[3], [0, 0, 0, 0])

        seqs = [DNASeq("ATACC"), DNASeq("AACC"), DNASeq("ATGAC")]
        m = p.distance_matrix(seqs)
        self.assertEqual(m[0], [0, 1, 2])
        self.assertEqual(m[1], [0, 0, 3])
        self.assertEqual(m[2], [0, 0, 0])

    def test_distance_matrix_aligned(self):
        p = PTree(sm, -1)
        seqs = [DNASeq("AT_ACC"), DNASeq("A__ACC"), DNASeq("ATGAC_")]
        m = p.distance_matrix_aligned(seqs)
        self.assertEqual(m[0], [0, 1, 2])
        self.assertEqual(m[1], [0, 0, 3])
        self.assertEqual(m[2], [0, 0, 0])

    def test_clustering(self):
        p = PTree(sm, -1)
        seqs = [DNASeq("A_CATATC_AT_"), DNASeq("A_GATATT_AG_"), DNASeq("AACAGATC_T__"), DNASeq("G_CAT__CGATT")]
        m = p.distance_matrix(seqs)
        c, t = p.clustering(m)
        self.assertEqual(c, "(((0,2),1),3)")
        self.assertIsNone(PTree.find_t(0, 1, []))
