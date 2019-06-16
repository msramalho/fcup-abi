import unittest
from bioseq import Blast
from bioseq import ProteinSeq
from bioseq.utils import *

db = "bioseq/resources/seqdump.txt"
fq = "bioseq/resources/source.fasta"
 
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
        self.assertEqual(b.forward_matching(q, q, (0,0), 0.5), len(q) - 3)
        self.assertEqual(b.forward_matching(ProteinSeq("PFMI"), ProteinSeq("PFMI"), (0,0), 0.5), 1)
        self.assertEqual(b.forward_matching(ProteinSeq("PFMI"), ProteinSeq("PFMI"), (3,3), 0.5), 0)
        self.assertEqual(b.extend_hit(ProteinSeq("PFMI"), ProteinSeq("PFMI"), (0,0), 0.5), 1)
        self.assertEqual(b.extend_hit(ProteinSeq("PFMI"), ProteinSeq("PFMI"), (3,3), 0.5), 1)

if __name__ == '__main__':
    unittest.main()
