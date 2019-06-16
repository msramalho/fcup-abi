from bioseq import BioSeq
from .matrix import Matrix
from .utils import GAP


class PTree:
    """Implements Phylogenetic Trees construction"""

    def __init__(self, sm, g):
        """Initialize the PTree class with the substitution matrix and gap penalty"""
        self.sm = sm
        self.g = g

    def distance_matrix(self, seqs):
        """Returns the distance matrix between all pairs in seqs"""
        m = Matrix(len(seqs))
        for i in range(len(seqs)):
            for j in range(i + 1, len(seqs)):
                score, traceback = seqs[i].global_align_multiple_solutions(seqs[j], self.sm, self.g)
                l, r = next(seqs[i].recover_global_align_multiple_solutions(seqs[j], traceback))
                m[i][j] = sum(c == GAP for c in l) + sum(c == GAP for c in r)
        return m.symmetric()
