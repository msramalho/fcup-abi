from collections import Counter


class BioSeq(object):
    """General purpose class for handling and manipulating biological sequences and perform operations on them."""
    valid_types = {"DNA", "RNA", "PROTEIN"}  #: What are the valid values for seq_type
    valid_tokens = {}  #: Which tokens a given biological sequence can have
    reverse = {}  #: rules for reverse_complement

    def __init__(self, sequence, seq_type):
        """Constructor for the Bio Sequence class"""
        self.sequence = sequence.upper()  #: The current sequence
        self.seq_type = seq_type.upper()  #: DNA, RNA, PROTEIN
        assert self._assert_seq_type(), "%s is not a valid sequence type(%s)" % (self.seq_type, BioSeq.valid_types)

    def frequency(self):
        """Calculates the relative frequency of each token in the sequence"""
        return {k: v / len(self) for k, v in Counter(self.sequence).items()}

    def gc_content(self):
        """Calculate the gc_content of the sequence"""
        return (self.sequence.count("G") + self.sequence.count("C")) / len(self)

    def reverse_complement(self):
        """Given a Biological sequence, generate the reverse complement sequence for the specific type. Uses the static variable reverse"""
        return "".join(type(self).reverse[x] for x in self.sequence)[::-1]

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

    def __repr__(self):
        """Class information visualization"""
        return self.__str__()

    def __getitem__(self, index):
        """get specific element of sequence"""
        return self.sequence[index]

    def __getslice__(self, start, stop, step=1):
        """get elements from [start to stop[ of sequence with a custom step"""
        return self.sequence[start:stop:step]
