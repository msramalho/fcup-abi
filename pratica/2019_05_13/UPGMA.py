from NumMatrix import NumMatrix
from HierarchicalClustering import HierarchicalClustering
from MySeq import MySeq
from PairwiseAlignment import PairwiseAlignment
from SubstMatrix import SubstMatrix

class UPGMA:

    def __init__(self, seqs, alseq):
        self.seqs = seqs
        self.alseq = alseq
        self.create_mat_dist()
        
    def create_mat_dist(self):
        # create distance matrix
        self.matdist = NumMatrix(#...)
        for i in range(len(self.seqs)):
            for j in range(i, len(self.seqs)):
                # retrieve the two sequences to align
                s1 = # ...
                s2 = # ...
                # align the sequences
                self.alseq.needleman_Wunsch(s1, s2)
                # recover the alignment
                alin = # ...
                ncd = 0
                # fill the matrix
                # count the number of different symbols in the alignment as the distance
                for k in range(len(alin)):
                    col = alin.column(k)
                    if (col[0] != col[1]): ncd += 1
                # set distance value in the matrix
                self.matdist.#...
                    
    def run(self):
        # create an object of the class HierarchicalClustering
        ch = # ...
        # execute the clustering algorithm
        t = # ...
        return t

def test():
    seq1 = MySeq("ATAGCGAT")    
    seq2 = MySeq("ATAGGCCT")    
    seq3 = MySeq("CTAGGCCC")
    seq4 = MySeq("CTAGGCCT")    
    sm = SubstMatrix()    
    sm.create_submat(1, -1, "ACGT")    
    alseq = PairwiseAlignment(sm, -2)    
    up  = UPGMA([seq1, seq2, seq3, seq4], alseq)    
    arv = up.run()    
    arv.print_tree() 

    

if __name__ == '__main__': 
    test()
    print()
    exercise1()
    
