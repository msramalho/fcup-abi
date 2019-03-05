import unittest
from bioseq import get_dict_of_dna_to_aminoacids, read_fasta


class TestUtils(unittest.TestCase):
    def test_get_dict_of_dna_to_aminoacids(self):
        self.assertDictEqual({'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'UGU': 'C', 'UGC': 'C', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'UUU': 'F', 'UUC': 'F', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'CAU': 'H', 'CAC': 'H', 'AUA': 'I', 'AUU': 'I', 'AUC': 'I', 'AAA': 'K', 'AAG': 'K', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'AUG': 'M', 'AAU': 'N', 'AAC': 'N', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'UGG': 'W', 'UAU': 'Y', 'UAC': 'Y', 'UAA': '_', 'UAG': '_', 'UGA': '_'}, get_dict_of_dna_to_aminoacids())

    def test_read_fasta(self):
        f = list(read_fasta("test/fasta_example.txt"))
        self.assertEqual(len(f), 16)
        self.assertEqual(str(f[0]), "('ORF number 1 in reading frame 1 on the direct strand extends from base 16 to base 165.', ATGATCAAGTTCATCAAAATCTTATTGGGGTGCCAAATAAACGTACCCTTGAATTTGCAAAATATTTGCAAAAACGTAATCAACATACCTGGATTCGTTATGTTGTGGTTCCTGGTTATACTGATAGCGATCACGATGTGCATTTATTAG)")

if __name__ == '__main__':
    unittest.main()
