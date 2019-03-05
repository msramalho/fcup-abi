import unittest
from bioseq import BioSeq


class TestBioSeqMethods(unittest.TestCase):
    def test_constructor(self):
        s = BioSeq("ATATAT", "DNA")
        self.assertRaises(Exception, BioSeq, "ATAT", "smth")

    def test_frequency(self):
        s = BioSeq("AACCCCCTGG", "DNA")
        self.assertDictEqual({"A": 20.0, "C": 50.0, "T": 10.0, "G": 20.0}, s.frequency())

    def test_assert_seq_type(self):
        s = BioSeq("ATATAT", "DNA")
        self.assertTrue(s._assert_seq_type())
        s.seq_type = "RNA"
        self.assertTrue(s._assert_seq_type())
        s.seq_type = "PROTEIN"
        self.assertTrue(s._assert_seq_type())
        s.seq_type = "ERROR"
        self.assertFalse(s._assert_seq_type())
    
    def test_str(self):
        s = BioSeq("ATATAT", "DNA")
        self.assertEqual("DNA: 'ATATAT'", str(s))

    def test_len(self):
        s = BioSeq("ATATAT", "DNA")
        self.assertEqual(6, len(s))


if __name__ == '__main__':
    unittest.main()
