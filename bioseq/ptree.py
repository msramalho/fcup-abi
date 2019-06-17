from bioseq import BioSeq
from .proteinseq import ProteinSeq
from .matrix import Matrix
from .utils import GAP
from .btree import BTree
from tqdm import tqdm


class PTree:
    """Implements Phylogenetic Trees construction"""

    def __init__(self, sm, g):
        """Initialize the PTree class with the substitution matrix and gap penalty"""
        self.sm = sm
        self.g = g
        self.ignore = 0

    def distance_matrix(self, seqs):
        """Returns the distance matrix between all pairs in seqs, assumes not aligned"""
        m = Matrix(len(seqs), value=self.ignore)
        with tqdm(total=sum(range(len(seqs)))) as pbar:
            for i in range(len(seqs)):
                for j in range(i + 1, len(seqs)):
                    score, traceback = seqs[i].global_align_multiple_solutions(seqs[j], self.sm, self.g)
                    l, r = next(seqs[i].recover_global_align_multiple_solutions(seqs[j], traceback))
                    m[i][j] = ProteinSeq(l).hamming_distance(ProteinSeq(r))
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

    def matrix_to_tuples(m):
        """Converts the superior triangular matrix into tuples of (row, col, distance) which is a simpler way of implementing UPGMA"""
        t = set()
        for i in range(len(m)):
            for j in range(i + 1, len(m[0])):
                t.add((i, j, m[i][j]))
        return t, set(range(len(m)))

    def clustering(self, m):
        trees = {i: BTree(i) for i in range(len(m))}
        t, s = PTree.matrix_to_tuples(m)
        while len(t) != 1:
            # get min
            mnt = min(t, key=lambda x: x[2])
            t -= {mnt}  # remove this tuple from t
            i, j, mn = mnt  # read the values
            trees["(%s,%s)" % (i,j)] = BTree.merge(trees[i], trees[j], mn)
            # remove i, j from s
            s -= {i, j}
            # calculate average for each remaining index
            for c in set(s):  # work on a copy of s, for each individual cluster
                x, y = PTree.find_t(i, c, t), PTree.find_t(j, c, t)  # find the relevant values
                t -= set([x, y])  # remove those tuples
                newt = PTree.merge_t(x, y, i, j, c)  # merge the tuples
                t.update([newt])  # add the new tuple
                s.add(newt[0])  # register the new cluster
        f = list(t)[0]
        return "(%s,%s)" % (f[0], f[1]), BTree.merge(trees[f[0]], trees[f[1]], mn)

    def find_t(i, c, t):
        # return [x for x in t if (x[0]==i and x[1]==c) or (x[0]==c and x[1]==i)]
        for x in t:
            if (x[0] == i and x[1] == c) or (x[0] == c and x[1] == i): return x
        return None

    def merge_t(x, y, i, j, c):
        return ("(%s,%s)" % (i, j), c, (x[2] + y[2]) / 2)
