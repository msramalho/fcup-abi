import unittest
from bioseq.utils import *


class TestUtils(unittest.TestCase):
    def test_get_dict_of_dna_to_aminoacids(self):
        self.assertDictEqual({'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'UGU': 'C', 'UGC': 'C', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'UUU': 'F', 'UUC': 'F', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'CAU': 'H', 'CAC': 'H', 'AUA': 'I', 'AUU': 'I', 'AUC': 'I', 'AAA': 'K', 'AAG': 'K', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'AUG': 'M', 'AAU': 'N', 'AAC': 'N', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'UGG': 'W', 'UAU': 'Y', 'UAC': 'Y', 'UAA': GAP, 'UAG': GAP, 'UGA': GAP}, get_dict_of_dna_to_aminoacids())

    def test_read_fasta(self):
        f = list(read_fasta("test/fasta_example.txt"))
        self.assertEqual(len(f), 4)
        self.assertEqual(str(f[0]), "('NP_004351.1 cadherin-1 isoform 1 preproprotein [Homo sapiens]', MGPWSRSLSALLLLLQVSSWLCQEPEPCHPGFDAESYTFTVPRRHLERGRVLGRVNFEDCTGRQRTAYFSLDTRFKVGTDGVITVKRPLRFHNPQIHFLVYAWDSTYRKFSTKVTLNTVGHHHRPPPHQASVSGIQAELLTFPNSSPGLRRQKRDWVIPPISCPENEKGPFPKNLVQIKSNKDKEGKVFYSITGQGADTPPVGVFIIERETGWLKVTEPLDRERIATYTLFSHAVSSNGNAVEDPMEILITVTDQNDNKPEFTQEVFKGSVMEGALPGTSVMEVTATDADDDVNTYNAAIAYTILSQDPELPDKNMFTINRNTGVISVVTTGLDRESFPTYTLVVQAADLQGEGLSTTATAVITVTDTNDNPPIFNPTTYKGQVPENEANVVITTLKVTDADAPNTPAWEAVYTILNDDGGQFVVTTNPVNNDGILKTAKGLDFEAKQQYILHVAVTNVVPFEAVSLTTSTATVTVDVLDVNEAPIFVPPEKRVEVSEDFGVGQEITSYTAQEPDTFMEQKITYRIWRDTANWLEINPDTGAISTRAELDREDFEHVKNSTYTALIIATDNGSPVATGTGTLLLILSDVNDNAPIPEPRTIFFCERNPKPQVINIIDADLPPNTSPFTAELTHGASANWTIQYNDPTQESIILKPKMALEVGDYKINLKLMDNQNKDQVTTLEVSVCDCEGAAGVCRKAQPVEAGLQIPAILGILGGILALLILILLLLLFLRRRAVVKEPLLPPEDDTRDNVYYYDEEGGGEEDQDFDLSQLHRGLDARPEVTRNDVAPTLMSVPRYLPRPANPDEIGNFIDENLKAADTDPTAPPYDSLLVFDYEGSGSEAASLSSLNSSESDKDQDYDYLNEWGNRFKKLADMYGGGEDD)")

    def test_module_exists(self):
        self.assertFalse(module_exists("adfkaodfn"))
        self.assertTrue(module_exists("collections"))

    def test_substitution_matrix(self):
        self.assertDictEqual({"AA": 2, "TT": 2, "AT": -3, "TA": -3}, substitution_matrix("AT", 2, -3))

    def test_read_substitution_matrix_file(self):
        sm = read_substitution_matrix_file("test/blosum62.mat")
        self.assertEqual(400, len(sm))
        for k in sm:
            self.assertEqual(sm[k], sm[k[::-1]])

    def test_score_column(self):
        sm = substitution_matrix("ATCG", 2, -3)
        self.assertEqual(2,  score_column("A", "A", sm, -5))
        self.assertEqual(-5, score_column("_", "A", sm, -5))
        self.assertEqual(-5, score_column("T", "_", sm, -5))
        self.assertEqual(-3, score_column("G", "C", sm, -5))

if __name__ == '__main__':
    unittest.main()
