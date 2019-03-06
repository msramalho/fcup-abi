from bioseq import *
import subprocess
import os


def sep(do_input=True):
    if do_input:
        input()
    print("-" * 60)

print("Hello!!\n")
print("This file will exemplify the capabilities of the developed project, this is an iterative process and you need to press enter to move from example to example, ok?")

sep()

print("You can create instances of BioSeq (with explicit sequence type) or one of the child classes DNASeq, RNASeq or ProteinSeq. \nEach has some common functions and some more specific. \nLet's create an instance of each...")
bs = BioSeq("ATCTCGTGCTGCTACG", "DNA")
print("Bioseq:", bs)
dna = DNASeq("ATCTGCTGCTC")
print("DNA:", dna)
rna = RNASeq("AUCGAUGCAUGCAC")
print("RNA:", rna)
pro = ProteinSeq("MNEPQLKHRASDYAQTQTQHYSEC_")
print("Protein:", pro)


sep()

print("The child classes validate the specific tokens in each sequence, so they will raise exceptions when user input is incorrect, they also take care of representing everything CAPS LOCK. Please now write a DNA sequence (valid or not) to check this behavior (example valid: 'ATCgtcCG', example invalid: 'ATCUGC')")

try:
    tmp = DNASeq(input())
    print("This sequence is valid", tmp)
except Exception as e:
    print("This sequence is invalid:", str(e))

sep()


print("Let us invoke some methods common to all classes because of BioSeq:")
print("Original sequence:", bs)
print("Frequency:", bs.frequency())
print("gc_content:", bs.gc_content())
print("pretty_print:")
bs.pretty_print()
filename = "temp.txt"
print("Save to file (%s), ok?" % filename)
input()
bs.save(filename)
print("Read from file (%s), ok? (you can go and change the sequence on the file, if you change the type it must also be valid)" % filename)
input()
bs.load(filename)
print("Read: ", bs)
os.remove(filename)
print("Deleted temporary file:", filename)

sep()
print("The special functions __str__, __repr__, __len__, __getitem__, __getslice__ have all been implemented so that we can access the real sequence directly from the class instances")


sep()
print("Let us invoke the specific DNA methods:")
print("Original DNA:", dna)
print("Transcript DNA (RNA):", dna.transcription())
print("translation:", dna.translation())
print("reverse_complement:", dna.reverse_complement())
print("codon_usage for aminoacid A:", dna.codon_usage("A"))
print("open_reading_frames:", list(dna.open_reading_frames()))

sep()
print("Let us invoke the specific RNA methods (most already seen for DNA, but it should be noted that in the DNA class the method involve a transcription and then they are the result of invoking on the RNA sequence, which is done here directly):")
print("Let us now use a longer RNA sequence so as to better see some ORFs:")
rna = RNASeq("AUGAAAUUAUGAAUGAGCCUCAGCUGAAGCAUCGCGCAUCAGACUACGCUCAGACUCAGACUCAGCAUUAUAGUGAAUGUUAAUAAAUAAAAUAA")
print("Original RNA:", rna)
print("translation:", rna.translation())
print("reverse_complement:", rna.reverse_complement())
print("codon_usage for aminoacid A:", rna.codon_usage("A"))
print("open_reading_frames:", list(rna.open_reading_frames()))


sep()

print("Now all the developed tests will be invoked, they thouroughly test the module and guarantee 100% coverage. ok?")
input()
subprocess.call("coverage run -m unittest discover -v".split())
