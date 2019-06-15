from bioseq import *

for f in read_fasta("bioseq/resources/source.fasta"):
    print(type(f[1]))