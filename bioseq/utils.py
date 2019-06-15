import bioseq
from pkgutil import iter_modules

GAP = "_"


def get_dict_of_dna_to_aminoacids():
    """Get a dictionary that maps RNA codons onto their respective proteins"""
    return {'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'UGU': 'C', 'UGC': 'C', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'UUU': 'F', 'UUC': 'F', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'CAU': 'H', 'CAC': 'H', 'AUA': 'I', 'AUU': 'I', 'AUC': 'I', 'AAA': 'K', 'AAG': 'K', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'AUG': 'M', 'AAU': 'N', 'AAC': 'N', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'UGG': 'W', 'UAU': 'Y', 'UAC': 'Y', 'UAA': GAP, 'UAG': GAP, 'UGA': GAP}


def read_fasta(filename):
    """Function that reads a FASTA file into a list of its ProteinSeq entries. Returns list of tuples (comment, ProteinSeq)"""
    with open(filename) as f:
        fastas = ("\n" + f.read().lstrip()).split("\n>")  # force split on descriptions
        for f in fastas[1:]:  # ignore the first empty
            p = f.split("\n")
            yield (p[0], bioseq.proteinseq.ProteinSeq("".join(p[1:])))


def module_exists(module_name):
    """Test if a module exists"""
    return module_name in (name for loader, name, ispkg in iter_modules())


def substitution_matrix(alphabet, match, mismatch):
    """Generate a substitution matrix from an alphabet, matcha and mismatch values"""
    return {x + y: match if x == y else mismatch for x in alphabet for y in alphabet}


def read_substitution_matrix_file(filename):
    """read substitution matrix from file """
    sm = {}
    with open(filename) as f:
        alphabet = next(f).split()
        for c in alphabet:
            for i, v in enumerate(map(int, next(f).split())):
                sm[c + alphabet[i]] = v
    return sm


def score_pos(c1, c2, sm, g):
    """score of a position (column). receives substituion_matrix and gap"""
    return g if c1 == GAP or c2 == GAP else sm[c1 + c2]


def kmer_generator(seq, k):
    """Given a sliceable object seq and a value for k, this function returns a generator for contiguous k-mers in seq"""
    for i in range(len(seq) - k + 1):
        yield seq[i:i + k]
