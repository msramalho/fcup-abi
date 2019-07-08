# globals for efficiency
DICT_AMINOACIDS = None


# return a valid dna sequence as an example
def dna_valid():
    return "ATGAAATTATGAATGAGCCTCAGCTGAAGCATCGCGCATCAGACTACGCTCAGACTCAGACTCAGCATTATAGTGAATGTTAATAAATAAAATAA"


# return an invalid dna sequence as an example
def dna_invalid():
    return "xGAAATTATGAATGAGCCTCAGCTGAAGCATCGCGCATCAGACTACGCTCAGACTCAGACTCAGCATTATAGTGAATGTTAATAAATAAAATAA"


# check if a dna sequence is valid
def is_valid(dna, codes=set(["A", "T", "C", "G"])):
    return all(x in codes for x in dna.upper())


# get a dict that translates codons of length 3 into aminoacids
def get_dict_of_aminoacids(filename="genetic_code.txt"):
    with open(filename) as f:
        return {p[0]: p[1] for p in [line.replace("\"", "").strip().split(" ") for line in f]}


# using get_dict_of_aminoacids translate a codon into an aminoacid with lazy loading
def translate_codon(cod):
    global DICT_AMINOACIDS
    DICT_AMINOACIDS = DICT_AMINOACIDS or get_dict_of_aminoacids()
    return DICT_AMINOACIDS[cod] if cod in DICT_AMINOACIDS else None


# read a dna sequence from a multiline file
def read_sequence(filename="example_Hinfluenzae.txt"):
    with open(filename) as f:
        return "".join(l.strip() for l in f)
