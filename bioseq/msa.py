from .utils import GAP
from collections import Counter
import copy
from tqdm import tqdm


class MSA:
    """Progressive implementation of Multiple Sequence Alignment"""

    def __init__(self, sm, g):
        """Initialize the MSA class with the substitution matrix and gap penalty"""
        self.sm = sm
        self.g = g

    def align(self, sequences):
        """Receive a list of sequences to align and return the MSA result"""
        seqs = [copy.deepcopy(s) for s in sequences]
        c = seqs[0]
        aligned = [c]
        klass = c.__class__
        with tqdm(total=len(seqs)-1) as pbar:
            for s in seqs[1:]:
                score, traceback = c.global_align_multiple_solutions(s, self.sm, self.g)
                c, s = next(c.recover_global_align_multiple_solutions(s, traceback))
                aligned = self.update_aligned_with_gaps(aligned, c)
                aligned.append(klass(s))  # add temp alignments to the list of processed
                c = self.consensus(aligned + [s], klass)
                pbar.update()
        return c, aligned

    def update_aligned_with_gaps(self, aligned, l):
        """include the new gaps in the previously aligned sequences so that the consensus function can use all the previous sequences and reduce error progagation"""
        gaps = [i for i in range(len(l)) if l[i] == GAP]
        for i, a in enumerate(aligned):
            for g in gaps: a.add_gap(g)
            aligned[i] = a
        return aligned

    def consensus(self, seqs, klass):
        """Calculate the consensus of all sequences with the new sequence"""
        res = ""
        for i in range(len(seqs[0])):
            mc = Counter(self.get_col(seqs, i)).most_common()
            mc = list(filter(lambda x: x[0]!=GAP, mc))
            mx = mc[0][1]
            mc = filter(lambda x: x[1]==mx, mc)
            mc = sorted(map(lambda x: x[0], mc))
            res += mc[0]
        return klass(res)

    def get_col(self, seqs, col):
        """return all the characters in column index of several sequences"""
        return [s[col] for s in seqs]
