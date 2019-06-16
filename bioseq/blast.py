from .utils import *
from tqdm import tqdm


class Blast:
    """Handle Blast searches"""

    def __init__(self, database, k):
        """create the database and do initial setup, k is the initial length to use for k-mers, database can be a list of ProteinSeq or a filename"""
        self.build_db(database, k)
        self.k = k

    def build_db(self, database, k):
        # check optional type for database
        if type(database) == str: database = list(read_fasta(database))
        # assert database characteristics
        assert type(database) == list, "Database argument must be a list"
        assert len(database), "At least one Protein sequence must be present"
        assert all(type(d) == bioseq.proteinseq.ProteinSeq for d in database), "All elements of database must be ProteinSeq"
        # build database kmers
        self.db = {seq: kmer_dict(seq, k) for seq in database}

    def search(self, query, top=10, threshold=0.5):
        """Perform a blast search of a given Protein sequence as query, return the best 'top' sequences"""
        # get kmer dict for query sequence
        # get hits
        hits = self.get_hits(query)
        # produce scores
        scores = self.get_hits_score(query, hits, threshold)
        # sort
        print(scores)
        # return first top
        pass

    def get_hits_score(self, query, hits, threshold):
        scores = defaultdict(int)
        with tqdm(total=sum(len(h) for h in hits.values())) as pbar:
            for seq, s_hits in hits.items():
                for h in s_hits:
                    scores[seq] = max(scores[seq], self.extend_hit(query, seq, h, threshold))
                    pbar.update()
        return scores

    def get_hits(self, query):
        """return a dict of {seq: list of tuples (query_index, seq_index)} where there are matches between database sequences and the query sequence"""
        hits = defaultdict(list)
        for seq, kd in self.db.items():
            if seq == query: continue
            for i, kmer in kmer_generator(query, self.k):
                if kmer in kd:
                    for db_index in kd[kmer]: hits[seq].append((i, db_index))
        return hits

    def extend_hit(self, query, seq, hit, threshold):
        """calculate the score of extending forwards and backwards, see forward_matching"""
        cf, forward = self.forward_matching(query, seq, hit, threshold)
        rhit = (len(query) - hit[0] - 1, len(seq) - hit[1] - 1)
        cb, backward = self.forward_matching(query[::-1], seq[::-1], rhit, threshold)
        return (cf + cb) / (forward + backward)

    def forward_matching(self, query, seq, hit, threshold):
        """advances on both sequences while ratio of matches is above threshold. returns the number of matches and the number of advances made"""
        qi, si = hit  # start of both sequences
        c = 0  # count matches while moving
        while qi + self.k < len(query) and si + self.k < len(seq):
            qi += 1; si += 1
            c += int(query[qi + self.k - 1] == seq[si + self.k - 1])
            # assert threshold condition
            if c / (qi - hit[0]) < threshold: break
        return c, qi - hit[0]
