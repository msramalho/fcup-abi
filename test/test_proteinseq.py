import unittest
from bioseq import ProteinSeq


class TestProteinSeq(unittest.TestCase):
    def test_constructor(self):
        s = ProteinSeq("MKL_MSLS_SIAHQTTLRLRLSIIVNVNK_N")
        self.assertEqual("MKL_MSLS_SIAHQTTLRLRLSIIVNVNK_N", s.sequence)
        self.assertEqual("PROTEIN", s.seq_type)
    
    def test_custom_assert_valid(self):
        s = ProteinSeq("MKL_MSLS_SIAHQTTLRLRLSIIVNVNK_N")
        self.assertTrue(s._assert_valid_sequence())
        s.sequence+="2"
        self.assertFalse(s._assert_valid_sequence())

    def test_gc_content(self):
        s = ProteinSeq("MKL_MSLS_SIAHQTTLRLRLSIIVNVNK_N")
        self.assertRaises(NotImplementedError, s.gc_content)
        
    def test_reverse_complement(self):
        s = ProteinSeq("MKL_MSLS_SIAHQTTLRLRLSIIVNVNK_N")
        self.assertRaises(NotImplementedError, s.reverse_complement)



if __name__ == '__main__':
    unittest.main()
