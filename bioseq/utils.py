import bioseq

def get_dict_of_dna_to_aminoacids():
    return {'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'UGU': 'C', 'UGC': 'C', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'UUU': 'F', 'UUC': 'F', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'CAU': 'H', 'CAC': 'H', 'AUA': 'I', 'AUU': 'I', 'AUC': 'I', 'AAA': 'K', 'AAG': 'K', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'AUG': 'M', 'AAU': 'N', 'AAC': 'N', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'UGG': 'W', 'UAU': 'Y', 'UAC': 'Y', 'UAA': '_', 'UAG': '_', 'UGA': '_'}


def read_fasta(filename):
    with open(filename) as f:
        fastas = ("\n" + f.read().lstrip()).split("\n>")  # force split on descriptions
        for f in fastas[1:]:  # ignore the first empty
            p = f.split("\n")
            yield (p[0], bioseq.proteinseq.ProteinSeq("".join(p[1:])))
