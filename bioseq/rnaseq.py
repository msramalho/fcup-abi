from bioseq import BioSeq


class RNASeq(BioSeq):
    """"Biological Sequence specific for RNA"""
    valid_tokens = {"A", "U", "C", "G"}

    def __init__(self, sequence):
        super().__init__(sequence, "RNA")
