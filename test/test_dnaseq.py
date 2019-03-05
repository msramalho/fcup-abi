import unittest
from bioseq import BioSeq, DNASeq


class TestBioSeqMethods(unittest.TestCase):
    def test_constructor(self):
        s = DNASeq("ATATAT")
        self.assertEqual("ATATAT", s.sequence)
        self.assertEqual("DNA", s.seq_type)
    
    def test_custom_assert_valid(self):
        s = DNASeq("ATATAT")
        self.assertTrue(s._assert_valid_sequence())
        s.sequence+="X"
        self.assertFalse(s._assert_valid_sequence())


if __name__ == '__main__':
    unittest.main()
