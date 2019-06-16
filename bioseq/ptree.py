from bioseq import BioSeq
from .matrix import Matrix
from .utils import GAP
from tqdm import tqdm

class PTree:
    """Implements Phylogenetic Trees construction"""

    def __init__(self, sm, g):
        """Initialize the PTree class with the substitution matrix and gap penalty"""
        self.sm = sm
        self.g = g
        self.ignore=0

    def distance_matrix(self, seqs):
        """Returns the distance matrix between all pairs in seqs, assumes not aligned"""
        m = Matrix(len(seqs), value=self.ignore)
        with tqdm(total=sum(range(len(seqs)))) as pbar:
            for i in range(len(seqs)):
                for j in range(i + 1, len(seqs)):
                    score, traceback = seqs[i].global_align_multiple_solutions(seqs[j], self.sm, self.g)
                    l, r = next(seqs[i].recover_global_align_multiple_solutions(seqs[j], traceback))
                    m[i][j] = l.hamming_distance(r)
                    pbar.update()
        return m

    def distance_matrix_aligned(self, seqs):
        """Returns the distance matrix between all pairs in seqs"""
        m = Matrix(len(seqs), value=self.ignore)
        with tqdm(total=sum(range(len(seqs)))) as pbar:
            for i in range(len(seqs)):
                for j in range(i + 1, len(seqs)):
                    m[i][j] = seqs[i].hamming_distance(seqs[j])
                    pbar.update()
        return m


    def clustering(self, m):
        mn = m.min(self.ignore)
