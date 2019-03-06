from collections import Counter
import re


class BioSeq(object):
    """General purpose class for handling and manipulating biological sequences and perform operations on them."""
    valid_types = {"DNA", "RNA", "PROTEIN"}  #: What are the valid values for seq_type
    valid_tokens = {}  #: Which tokens a given biological sequence can have
    reverse = {}  #: rules for reverse_complement

    def __init__(self, sequence, seq_type):
        """Constructor for the Bio Sequence class"""
        self.sequence = sequence.upper()  #: The current sequence
        self.seq_type = seq_type.upper()  #: DNA, RNA, PROTEIN
        self._assert_seq_type()

    def frequency(self):
        """Calculates the relative frequency of each token in the sequence"""
        return {k: v / len(self) for k, v in Counter(self.sequence).items()}

    def gc_content(self):
        """Calculate the gc_content of the sequence"""
        return (self.sequence.count("G") + self.sequence.count("C")) / len(self)

    def reverse_complement(self):
        """Given a Biological sequence, generate the reverse complement sequence for the specific type. Uses the static variable reverse"""
        return "".join(type(self).reverse[x] for x in self.sequence)[::-1]

    def pretty_print(self):
        """Outputs a prettified and informative string onto the console, describing the current sequence"""
        print("""Sequence Type:   %s\nLength:          %d\nToken frequency: %s\nSequence:\n%s""" % (
            self.seq_type, len(self), self.frequency(), "\n".join(self.sequence[i:i + 80] for i in range(0, len(self), 80))))

    def save(self, filename):
        """Save an object to a file"""
        with open(filename, "w") as fout:
            fout.write("%s\n%s" % (self.seq_type, self.sequence))

    def load(self, filename):
        """Load object from file"""
        with open(filename) as fin:
            self.seq_type = fin.readline().strip()
            self.sequence = fin.readline().strip()

    def _assert_valid_sequence_regex(self):
        """assert that all the tokens in the sequence are valid for that sequence type using a regex operator"""
        tokens = "".join(type(self).valid_tokens)
        pattern = "^[%s]*$" % (tokens + tokens.lower())
        assert bool(re.search(pattern, self.sequence)), "%s is not a valid sequence(%s)" % (self.sequence, type(self).valid_tokens)

    def _assert_valid_sequence(self):
        """assert that all the tokens in the sequence are valid for that sequence type"""
        assert all(x in type(self).valid_tokens for x in self.sequence), "%s is not a valid sequence(%s)" % (self.sequence, type(self).valid_tokens)

    def _assert_seq_type(self):
        """Assert a proper valid has been chosen for seq_type"""
        assert self.seq_type in BioSeq.valid_types, "%s is not a valid sequence type(%s)" % (self.seq_type, BioSeq.valid_types)

    def __len__(self):
        """Get the length of the sequence - number of tokens"""
        return len(self.sequence)

    def __str__(self):
        """Class information visualization"""
        return self.sequence

    def __repr__(self):
        """Class information visualization"""
        return self.__str__()

    def __getitem__(self, index):
        """get specific element of sequence"""
        return self.sequence[index]

    def __getslice__(self, start, stop, step=1):
        """get elements from [start to stop[ of sequence with a custom step"""
        return self.sequence[start:stop:step]
