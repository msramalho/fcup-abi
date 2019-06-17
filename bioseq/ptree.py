from bioseq import BioSeq
from .proteinseq import ProteinSeq
from .matrix import Matrix
from .utils import GAP
from .btree import BTree
from tqdm import tqdm
from Bio import Phylo
import uuid
import os, io

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
        self.btree = BTree.merge(trees[f[0]], trees[f[1]], mn)
        self.tree = self.phylo_tree()
        return "(%s,%s)" % (f[0], f[1]), self.btree
    
    def phylo_tree(self):
        tmp = "%s.xml" % uuid.uuid4().hex
        with open(tmp, "w") as o: o.write(self.btree.to_xml())
        t = Phylo.read(tmp,'phyloxml')
        os.remove(tmp)
        return t

    def find_t(i, c, t):
        """find the tuple with i and c as cluster ids, in whatever order"""
        for x in t:
            if (x[0] == i and x[1] == c) or (x[0] == c and x[1] == i): return x

    def merge_t(x, y, i, j, c):
        """return the merged cluster from two instances, averaging the distances"""
        return ("(%s,%s)" % (i, j), c, (x[2] + y[2]) / 2)

    def draw(self): # pragma: no cover
        if "tree" not in dir(self): return "call clustering first"
        Phylo.draw(self.tree)

    def __repr__(self): # pragma: no cover
        return self.__str__()

    def __str__(self): # pragma: no cover
        if "tree" not in dir(self): return "call clustering first"
        with io.StringIO() as output:
            Phylo.draw_ascii(self.tree, file=output)
            return output.getvalue()

