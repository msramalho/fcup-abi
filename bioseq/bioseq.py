import re
from operator import itemgetter
from collections import Counter
from .matrix import Matrix
from .utils import *

# traceback settings
TERMINATION = -1
DIAGONAL = 0
VERTICAL = 1
HORIZONTAL = 2


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

    def hamming_distance(self, seq):
        """implements a hamming distance calculator"""
        assert len(self) == len(seq), "Hamming distance is only possible between same length sequences"
        return sum(self[i] != seq[i] for i in range(len(self)))

    def dot_plot(self, seq):
        """create a dotplot between two sequences. Returns a Matrix"""
        m = Matrix(len(self), len(seq))
        m.apply(lambda _, i, j: int(self[i] == seq[j]))
        return m

    def display_dot_plot(self, m):  # pragma: no cover
        """Display the dot_plot in a plot if matplotlib is installed"""
        m.display()

    def score_seq(self, seq, sm, g):
        """Calculate the score of aligning two sequences, given a gap and substitution matrix"""
        assert len(self) == len(seq), "Sequences should have the same length"
        return sum(score_pos(self[i], seq[i], sm, g) for i in range(len(self)))

    def score_affine_gap(self, seq, sm, g, r):
        """Calculate the affine gap score, using the gap penalty and the keep gap penalty"""
        assert len(self) == len(seq), "Sequences should have the same length"
        gap = False
        res = 0
        for i in range(len(self)):
            if self[i] == GAP or seq[i] == GAP:
                if not gap: gap = True; r += g
                else: res += r
            else:
                gap = False
                res += sm[self[i] + seq[i]]
        return res

    def global_align_multiple_solutions(self, seq, sm, g):
        """Needleman–Wunsch"""
        # create the score and traceback matrices
        s = Matrix(len(self) + 1, len(seq) + 1)
        t = Matrix(len(self) + 1, len(seq) + 1, [TERMINATION])
        # set the row and col to gaps
        s[0] = [i * g for i in range(len(seq) + 1)]
        s.set_col(0, [i * g for i in range(len(self) + 1)])
        # set the row and col to default directions
        t[0] = [[HORIZONTAL]] * (len(seq) + 1)
        t.set_col(0, [[VERTICAL]] * (len(self) + 1))
        # fill score matrix
        for i in range(1, len(self) + 1):  # for each row
            for j in range(1, len(seq) + 1):  # for each col
                d = [s[i - 1][j - 1] + sm[self[i - 1] + seq[j - 1]], s[i - 1][j] + g, s[i][j - 1] + g]
                s[i][j] = max(d)  # set the score
                t[i][j] = [k for k, v in enumerate(d) if v == s[i][j]]  # set directions
        return s, t

    def _recover_global_dfs(self, seq, t, i, j):
        """helper function that uses DFS to build multiple results"""
        if i == 0 and j == 0: return [("","")]
        chains = []
        res = ["", ""]
        for step in t[i][j]:
            if step == HORIZONTAL: j-=1; res[0] += GAP; res[1] += seq[j]
            elif step == VERTICAL: i-=1; res[1] += GAP; res[0] += self[i]
            else: i-=1; j-=1; res[0] += self[i]; res[1] += seq[j]
            for a, b in self._recover_global_dfs(seq, t, i, j):
                chains.append((res[0] + a, res[1] + b))
        return chains

    def recover_global_align_multiple_solutions(self, seq, t):
        """Given two sequences and their Score and Traceback global alignment functions, return all the solutions"""
        return [(a[::-1], b[::-1]) for a, b in self._recover_global_dfs(seq, t, len(self), len(seq))]

    def local_align_multiple_solutions(self):
        """Smith–Waterman"""
        pass

    def recover_local_align_multiple_solutions(self):
        """"""
        pass

    def compare_pairwise_global_align(self):
        """"""
        pass

    def compare_pairwise_local_align(self):
        """"""
        pass

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
