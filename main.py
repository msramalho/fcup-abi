from bioseq import *

sm = read_substitution_matrix_file("test/blosum62.mat")



s1 = BioSeq("GATTACA", "PROTEIN")
s2 = BioSeq("GCATGCT", "PROTEIN")
sm = substitution_matrix("ATCG", 1, -1)

# print(sm)
s, t = s1.global_align_multiple_solutions(s2, sm, -1)
print(s)
print(t)
for pair in s1.recover_global_align_multiple_solutions(s2, t):
    print(pair)




sm = read_substitution_matrix_file("test/blosum62.mat")
s1 = BioSeq("HGWAG", "PROTEIN")
s2 = BioSeq("PHSWG", "PROTEIN")
s, t = s1.local_align_multiple_solutions(s2, sm, -8)
print(s)
print(t)
for pair in s1.recover_local_align_multiple_solutions(s2, t, s):
    print(pair)
# print(s)
# print(t)
# print(s2.recover_global_align_multiple_solutions(s1, t))

sm = read_substitution_matrix_file("test/blosum62.mat")
from_fasta = {x[0]:x[1] for x in read_fasta("test/protein_sequences.fas")}
s1 = from_fasta["sp|C1F111"]
s2 = from_fasta["sp|B7JC18"]
print(len(s1), len(s2))
s, t = s1.local_align_multiple_solutions(s2, sm, -8)
print(s.max())
for pair in s1.recover_local_align_multiple_solutions(s2, t, s):
    print(pair)

print(BioSeq.compare_pairwise_global_align(from_fasta.values(), sm, -3, list(from_fasta.keys())))
print(BioSeq.compare_pairwise_local_align(from_fasta.values(), sm, -3, list(from_fasta.keys())))