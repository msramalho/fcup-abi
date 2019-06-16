from bioseq import BioSeq
from .utils import GAP


class ProteinSeq(BioSeq):
    """Biological Sequence specific for PROTEINS"""
    valid_tokens = {'P', 'F', 'M', 'I', 'R', 'K', 'A', 'L', 'V', 'Q', 'E', 'C', 'N', 'W', 'H', 'S', 'T', 'D', 'G', 'Y', GAP}

    def __init__(self, sequence, name=None):
        super().__init__(sequence, "PROTEIN")
        self._assert_valid_sequence()
        if name: self.name = name

    def gc_content(self):
        raise NotImplementedError("The class ProteinSeq does not have a meaningful implementation for gc_content")

    def reverse_complement(self):
        raise NotImplementedError("The class ProteinSeq does not have a meaningful implementation for reverse_complement")
