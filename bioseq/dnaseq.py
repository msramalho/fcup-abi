from bioseq import BioSeq
from .rnaseq import RNASeq
from .utils import GAP


class DNASeq(BioSeq):
    """Biological Sequence specific for DNA"""
    valid_tokens = {"A", "T", "C", "G", GAP}
    reverse = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def __init__(self, sequence):
        super().__init__(sequence, "DNA")
        self._assert_valid_sequence()

    def transcription(self):
        """Performs a transcription of the original sequence from DNA to RNA and returns a RNASeq object"""
        return RNASeq(self.sequence.replace("T", "U"))

    def translation(self, start=0):
        """Given a DNA sequence, transcript it to RNA and call the respective translation method for RNA sequences. Starts codon split at position **start** (default is 0). Returns a ProteinSeq instance. """
        return self.transcription().translation(start)

    def codon_usage(self, aminoacid):
        """Given a sequence,  transcript it to RNA and call the respective codon_usage to get the frequency of each codon in it that maps to a specific aminoacid. Returns a dict. """
        return self.transcription().codon_usage(aminoacid)

    def open_reading_frames(self):
        """Given a sequence, get all the possible open reading frames for each reading frame after transcripting to RNASeq"""
        for orf in self.transcription().open_reading_frames():
            yield orf
