from bioseq import *
import subprocess
import os


def sep(do_input=True):
    if do_input:
        input()
    print("-" * 60)


print("Hello!!\n")
print("This file will exemplify the capabilities of the developed project, this is an iterative process and you need to press enter to move from example to example [according to the project script requirements], ok?")

sep()

print("This script shows the new functionality concerning protein alignment, the first part can be seen with 'python p1_run_me.py'")

sep()

print("Section A: global_align_multiple_solutions")
s1 = ProteinSeq("PHSWG")
s2 = ProteinSeq("HGWAG")
print("using sequences:", s1, s2)
input()
sm = read_substitution_matrix_file("test/blosum62.mat")
print("using the BLOSUM62 matrix", sm)
input()
s, t = s1.global_align_multiple_solutions(s2, sm, -8)
print("the score and traceback of global_align_multiple_solutions with gap=-8 (example from the slides) is")
print(s)
print(t)
print("OPTIMAL alignment: ", s.last())
input()

sep()

print("Section A: recover_global_align_multiple_solutions for the previous alignment:")
for x in s1.recover_global_align_multiple_solutions(s2, t):
    print(x)

sep()
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

sep()

print("Section C1: protein_sequence.fas BLOSUM62 gap=-3")
from_fasta = {x[0]:x[1] for x in read_fasta("test/protein_sequences.fas")}
print("Reading from 'protein_sequence.fas' yields the following dict NAME->BioSeq", from_fasta)
input()
print("Calculating at least a pair with multiple alignments")
found = False
for k, v in from_fasta.items():
    if found: break
    for k1, v1 in from_fasta.items():
        print("Testing ", k, k1)
        if k == k1: continue
        s, t = v.global_align_multiple_solutions(v1, sm, -3)
        rec = list(v.recover_global_align_multiple_solutions(v1, t))
        if len(rec) > 1:
            print("FOUND with %d optimal alignments" % len(rec))
            input()
            print(rec)
            found = True
            break
input()


print("Section C2:  GATTACA and GCATGCT with match = 1, mismatch = -1 and gap = -1")

s1 = DNASeq("GATTACA")
s2 = DNASeq("GCATGCT")
print("using sequences:", s1, s2)
input()
sm = substitution_matrix("ACTG", 1, -1)
print("substitution matrix: ", sm)
input()
print("global alignment matrices:")
s, t = s1.global_align_multiple_solutions(s2, sm, -1)
print(s)
print(t)
input()
rec = list(s1.recover_global_align_multiple_solutions(s2, t))
print(len(rec), "alignments:", rec)

sep()

print("Section D: compare_pairwise_global_align")
sm = read_substitution_matrix_file("test/blosum62.mat")
print("Using the sequences from 'protein_sequences'")
print(BioSeq.compare_pairwise_global_align(from_fasta.values(), sm, -3, list(from_fasta.keys())))
input()
print("Section D: compare_pairwise_local_align")
print(BioSeq.compare_pairwise_local_align(from_fasta.values(), sm, -3, list(from_fasta.keys())))

print("As you saw, this takes a while so it was optimized in the sense that only the triangular matrice values are calculated and the remaining values are filled using the apply method of the Matrix class :)")

sep()

print("Some additional functionality was also developed!")
print("The Matrix class that has been used throughout the previous demonstration.")
print("Let's consider the sequences:")
s1 = DNASeq("ATAT")
s2 = DNASeq("AGCT")
print(s1, s2)
input()

print("Dotplot for a matrix")
dp = s1.dot_plot(s2)
print(dp)
input()


print("Dotplots can also be visualized (only works if matplotlib is installed)")
dp.display()
input()

print("Dotplots can also be saved to a file, let's send it to 'temp.png'")
dp.display("temp.png")
input()

print("Hamming distance between sequences")
print("hamming distance:", s1, s2, s1.hamming_distance(s2))

sep()

print("Now all the developed tests will be invoked, they thouroughly test the module and guarantee 100% coverage. ok?")
input()
subprocess.call("coverage run -m unittest discover -v".split())
