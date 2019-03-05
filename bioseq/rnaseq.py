from bioseq import BioSeq
from .utils import get_dict_of_dna_to_aminoacids
from .proteinseq import ProteinSeq


class RNASeq(BioSeq):
    """Biological Sequence specific for RNA"""
    valid_tokens = {"A", "U", "C", "G"}
    reverse = {"A": "U", "U": "A", "C": "G", "G": "C"}
    dict_of_aminoacids = get_dict_of_dna_to_aminoacids()

    def __init__(self, sequence):
        super().__init__(sequence, "RNA")

    def translation(self, start=0):
        """Given an RNA sequence, translate its 3-token codons into the respective aminoacids. Starts codon split at position **start** (default is 0). Returns a ProteinSeq instance."""
        return ProteinSeq("".join(type(self).dict_of_aminoacids[self.sequence[i:i + 3]] for i in range(start, len(self) - 2, 3)))

    def codon_usage(self, aminoacid):
        """Given a sequence, calculate the frequency of each codon in it that maps to a specific aminoacid. Returns a dict. """
        codons = [k for k, v in type(self).dict_of_aminoacids.items() if v == aminoacid]  # get relevant codons
        rna_codons = [self.sequence[i:i + 3] for i in range(0, len(self) - 2, 3)]  # get all the codons in the dna sequence
        return {codon: rna_codons.count(codon) / len(rna_codons) for codon in codons}  # output the frequency
