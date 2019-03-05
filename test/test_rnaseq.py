import unittest
from bioseq import RNASeq


class TestRNASeqMethods(unittest.TestCase):
    def test_constructor(self):
        s = RNASeq("AUGCGAU")
        self.assertEqual("AUGCGAU", s.sequence)
        self.assertEqual("RNA", s.seq_type)
    
    def test_custom_assert_valid(self):
        s = RNASeq("AUGCGAU")
        self.assertTrue(s._assert_valid_sequence())
        s.sequence+="T"
        self.assertFalse(s._assert_valid_sequence())

    def test_reverse_complement(self):
        s = RNASeq("AUCG")
        self.assertEqual("CGAU", s.reverse_complement())

if __name__ == '__main__':
    unittest.main()
