from bioseq import BioSeq


class DNASeq(BioSeq):
    """"Biological Sequence specific for DNA"""
    valid_tokens = {"A", "T", "C", "G"}
    reverse = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def __init__(self, sequence):
        super().__init__(sequence, "DNA")
