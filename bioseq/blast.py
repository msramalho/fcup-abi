from .utils import *


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
        print(type(database[0]))
        assert all(type(d) == bioseq.proteinseq.ProteinSeq for d in database), "All elements of database must be ProteinSeq"
        # build database kmers
        self.db = {seq: kmer_dict(seq, k) for seq in database}

    def search(self, query, top=10):
        """Perform a blast search of a given Protein sequence as query, return the best 'top' sequences"""
        # get kmer dict for query sequence
        # get hits
        hits = self.get_hits(query)
        print("hits:", hits)
        # extend hits
        # produce score
        # sort
        # return first top
        pass

    def get_hits(self, query):
        hits = defaultdict(list)
        for seq, kd in self.db.items():
            print(".....")
            if seq == query: continue
            for i, kmer in kmer_generator(query, self.k):
                if kmer in kd: 
                    for db_index in kd[kmer]: hits[seq].append((i, db_index))
            break
        return hits
