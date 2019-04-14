from bioseq import *





print("Section B: local_align_multiple_solutions")

s1 = ProteinSeq("PHSWG")
s2 = ProteinSeq("HGWAG")
print("using sequences:", s1, s2)
input()
sm = read_substitution_matrix_file("test/blosum62.mat")
print("using the BLOSUM62 matrix, as before")
input()
s, t = s1.local_align_multiple_solutions(s2, sm, -8)
print("the score and traceback of global_align_multiple_solutions with gap=-8 (example from the slides) is")
print(s)
print(t)
print("OPTIMAL alignment: ", s.last())
input()

sep()

print("Section B: recover_local_align_multiple_solutions")
for x in s1.recover_local_align_multiple_solutions(s2, t, s):
    print(x)
