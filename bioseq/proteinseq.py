from bioseq import BioSeq


class ProteinSeq(BioSeq):
    """"Biological Sequence specific for PROTEINS"""
    valid_tokens = {'P', 'F', '_', 'M', 'I', 'R', 'K', 'A', 'L', 'V', 'Q', 'E', 'C', 'N', 'W', 'H', 'S', 'T', 'D', 'G', 'Y'}

    def __init__(self, sequence):
        super().__init__(sequence, "PROTEIN")
