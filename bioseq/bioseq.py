from collections import Counter


class BioSeq(object):
    """The summary line for a class docstring should fit on one line."""
    valid_types = {"DNA", "RNA", "PROTEIN"}
    valid_tokens = {}
    # valid_ami = {'P', 'F', '_', 'M', 'I', 'R', 'K', 'A', 'L', 'V', 'Q', 'E', 'C', 'N', 'W', 'H', 'S', 'T', 'D', 'G', 'Y'}
    # valid_rna = {"A", "U", "C", "G"}

    def __init__(self, sequence, seq_type):
        """Constructor for the Bio Sequence class"""
        self.sequence = sequence  #: The current sequence
        self.seq_type = seq_type.upper()  #: DNA, RNA, PROTEIN
        assert self._assert_seq_type(), "%s is not a valid sequence type(%s)" % (self.seq_type, BioSeq.valid_types)

    def frequency(self):
        """Calculates the relative frequency of each token in the sequence"""
        return {k: 100 * v / len(self) for k, v in Counter(self.sequence).items()}

    def _assert_valid_sequence(self):
        """assert that all the tokens in the sequence are valid for that sequence type"""
        return all(x in type(self).valid_tokens for x in self.sequence)

    def _assert_seq_type(self):
        """Assert a proper valid has been chosen for seq_type"""
        return self.seq_type in BioSeq.valid_types

    def __len__(self):
        """Get the length of the sequence - number of tokens"""
        return len(self.sequence)

    def __str__(self):
        """Class information visualization"""
        return "%s: '%s'" % (self.seq_type, self.sequence)
