from .utils import GAP


class MSA:
    """Progressive implementation of Multiple Sequence Alignment"""

    def __init__(self, sm, g):
        """Initialize the MSA class with the substitution matrix and gap penalty"""
        self.sm = sm
        self.g = g

    def align(self, seqs):
        """Receive a list of sequences to align and return the MSA result"""
        c = seqs[0]
        for s in seqs[1:]:
            score, traceback = c.global_align_multiple_solutions(s, self.sm, self.g)
            l, r = next(c.recover_global_align_multiple_solutions(s, traceback))
            c = self.consensus(l, r)
        return c

    def consensus(self, l, r):
        """Calculate the consensus of two sequences"""
        res = ""
        for i in range(len(l)):
            res += self.consensus_col(l[i], r[i])
        return l.__class__(res)

    def consensus_col(self, c1, c2):
        """Given two columns, return the most frequent, excluding gaps or the lexicographically smaller"""
        if c1 == GAP: return c2
        if c2 == GAP: return c1
        if c1 == c2: return c1
        if c1 < c2: return c1
        else: return c2
