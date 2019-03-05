from bioseq import BioSeq
from .rnaseq import RNASeq


class DNASeq(BioSeq):
    """Biological Sequence specific for DNA"""
    valid_tokens = {"A", "T", "C", "G"}
    reverse = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def __init__(self, sequence):
        super().__init__(sequence, "DNA")

    def transcription(self):
        """Performs a transcription of the original sequence from DNA to RNA and returns a RNASeq object"""
        return RNASeq(self.sequence.replace("T", "U"))

    def translation(self, start=0):
        """Given a DNA sequence, transcript it to RNA and call the respective translation method for RNA sequences. Starts codon split at position **start** (default is 0). Returns a ProteinSeq instance. """
        return self.transcription().translation(start)

    def codon_usage(self, aminoacid):
        """Given a sequence,  transcript it to RNA and call the respective codon_usage to get the frequency of each codon in it that maps to a specific aminoacid. Returns a dict. """
        return self.transcription().codon_usage(aminoacid)
