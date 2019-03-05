from bioseq import BioSeq


class DNASeq(BioSeq):
    """"Biological Sequence specific for DNA"""
    valid_tokens = {"A", "T", "C", "G"}

    def __init__(self, sequence):
        super().__init__(sequence, "DNA")
