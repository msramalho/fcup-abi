import io
import unittest
import unittest.mock
from bioseq import BioSeq


class TestBioSeq(unittest.TestCase):
    def test_constructor(self):
        self.maxDiff = None
        s = BioSeq("ATATat", "DNa")
        self.assertIsInstance(s, BioSeq)
        self.assertEqual("ATATAT", s.sequence)
        self.assertEqual("DNA", s.seq_type)
        self.assertRaises(Exception, BioSeq, "ATAT", "smth")

    def test_frequency(self):
        s = BioSeq("AACCCCCTGG", "DNA")
        self.assertDictEqual({"A": 0.2, "C": 0.5, "T": 0.1, "G": 0.2}, s.frequency())

    def test_gc_content(self):
        s = BioSeq("AACCCCCTGG", "DNA")
        self.assertEqual(0.70, s.gc_content())
        s.sequence = "AA"
        self.assertEqual(0.0, s.gc_content())
        s.sequence = "GC"
        self.assertEqual(1.0, s.gc_content())

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
        self.assertEqual("ATATAT", str(s))

    def test_repr(self):
        s = [BioSeq("ATATAT", "DNA")]
        self.assertEqual("[ATATAT]", str(s))

    def test_getitem(self):
        s = BioSeq("ATATAT", "DNA")
        self.assertEqual("A", s[0])
        self.assertEqual("T", s[-1])

    def test_getslice(self):
        s = BioSeq("ATATAT", "DNA")
        self.assertEqual("A", s[0:1])
        self.assertEqual("T", s[-1:])
        self.assertEqual("AT", s[0:2])
        self.assertEqual("A", s[0:1:1])
        self.assertEqual("T", s[-1::1])
        self.assertEqual("AAA", s[0::2])
        # explicit invocation of function
        self.assertEqual("A", s.__getslice__(0, 1))
        self.assertEqual("T", s.__getslice__(-1, len(s)))
        self.assertEqual("AT", s.__getslice__(0, 2))
        self.assertEqual("A", s.__getslice__(0, 1, 1))
        self.assertEqual("T", s.__getslice__(-1, len(s), 1))
        self.assertEqual("AAA", s.__getslice__(0, len(s), 2))

    def test_len(self):
        s = BioSeq("ATATAT", "DNA")
        self.assertEqual(6, len(s))

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def _test_pretty_print(self, bioseq, mock_stdout):
        bioseq.pretty_print()
        output = mock_stdout.getvalue()
        self.assertEqual(len(output), 186)
        self.assertIn(bioseq.sequence, output)

    def test_pretty_print(self):
        s = BioSeq("ATATACAGATGAT", "DNA")
        self._test_pretty_print(s)


if __name__ == '__main__':
    unittest.main()
