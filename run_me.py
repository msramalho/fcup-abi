#!/usr/bin/env python
# coding: utf-8

import subprocess
import os


def sep(do_input=True):
    print("-" * 60)
    if do_input: input()

# In[1]:


from bioseq import *

print("This run_me.py executes all the functions present in the Jupyter Notebook in the console. For a clearer display run the notebook OR check the run_me.html which is a static file with the execution result from the Jupyter Notebook. ok?")
input()
print("Furthermore, to see all the visualizations of trees and graphs some libraries are required, these are specified in the 'requirements.txt' file and the use of a virtual environment is advised when installing them. To install simply do 'pip install -r requirements.txt'.\n\nAfter this you are good to go!")
sep()

print("Load the substitution matrices for further use")

# In[9]:


blosum = read_substitution_matrix_file("bioseq/resources/blosum62.mat")
sm = substitution_matrix("ATGC_", 1, -1)
print("blosum", blosum)
input()
print("ATGC_ substitution matrix:", sm)

sep()


print("BLAST SECTION")
print("Perform blast on the provided database and sequence")
input()

# In[4]:


print("initialize (will load the database with k=3)")
b = Blast("bioseq/resources/seqdump.txt", 3)


# In[5]:


print("Items in the db: ", len(b.db))
print("Value to use for k in the kmers:", b.k)
input()

# In[6]:


print('Load the query sequence from the "source.fasta" file')
q = next(read_fasta("bioseq/resources/source.fasta"))
print(q)
input()


# In[7]:


print("Perform BLAST search over the previously loaded database")
print("By default returns top=10 with a search treshold=0.5")
similar = b.search(q)
input()


# In[8]:


print("Sorted scores:")
for seq, score in similar:
    print(seq.name, ":", score)
input()

# In[9]:


print("get a list with all the sequences (including the query sequence) for the next step: MSA")
blast_seqs = [x[0] for x in similar] + [q]
print(blast_seqs)
sep()

print("Multiple Sequence Alignment Section")

# In[11]:


print("create an MSA instance with blosum62 and a gap penalty of -1")
m = MSA(blosum, -1)
input()


# In[12]:

print("perform the MSA on the blast sequences")
consensus, aligned_seqs = m.align(blast_seqs)
input()


# In[13]:


print("Consensus:", consensus)
input()

# In[14]:


print("Aligned sequences:", aligned_seqs)
input()


print("MSA with a simpler (easier to check visually) example")

# In[15]:


print("Using Blosum62")
seqs = [ProteinSeq("PHWAS"), ProteinSeq("HWASW"), ProteinSeq("HPHWA")]
print("Original sequences are:", seqs)
input()
print("align result (consensus, aligned seqs):", m.align(seqs))
input()

# In[16]:


print("Using simpler substitution matrix and with the example from the practical lessons")
ms = MSA(sm, -1)


# In[17]:

seqs = [DNASeq("ATAGC"), DNASeq("AACC"), DNASeq("ATGAC")]
print("Original sequences used in class are:", seqs)
input()
print("align result (consensus, aligned seqs):", ms.align(seqs))
sep()

print("Phylogenetic Trees")
input()
print("1st: Using a simple example (from the classes)")

# In[18]:


print("simple sm and gap penalty of -1")
ps = PTree(sm, -1)
input()

# In[19]:


print("first produce the distance matrix")
seqs = [DNASeq("A_CATATC_AT_"), DNASeq("A_GATATT_AG_"), DNASeq("AACAGATC_T__"), DNASeq("G_CAT__CGATT")]
print("sequences are:", seqs)
input()
mxs = ps.distance_matrix(seqs)
print(mxs)
input()

# In[20]:


print("perform UPGMA clustering")
c, t = ps.clustering(mxs)
print("simple clustering representation", c)
input()


# In[21]:


print("ASCII tree representation")
print(ps)
input()


# In[22]:


print("drawing the tree with Phylo:")
ps.draw()
input()


# In[23]:


print("You can also check Phylo's tree text view with ps.tree")
print(ps.tree)
input()

sep()

print("2nd: Using the BLAST sequences (takes longer to execute ~200s)")
print("This takes longer to execute, but you can use the progress bar to keep track of time")
input()

# In[24]:


print("blosum and gap penalty of -1")
p = PTree(blosum, -1)
input()

# In[25]:

print("get distance matrix")
mx = p.distance_matrix(blast_seqs)
print(mx)
input()


# In[26]:
# This piece of code is only used to extract the contents within brackets from the species names
# so as to get more readable visualizations
from copy import deepcopy
def get_seq_with_name(x):
    y = deepcopy(x)
    n = y.name
    y.name = n[n.index("[")+1:n.index("]")]
    return y
blast_seqs_names = list(map(get_seq_with_name, blast_seqs))

print("perform UPGMA clustering")
c, t = p.clustering(mx, blast_seqs_names)
print("simple clustering representation", c)
input()


# In[27]:


print("ASCII tree representation")
print(p)
input()


# In[28]:


print("drawing the tree with Phylo:")
p.draw()
input()


# In[29]:


print("You can also check Phylo's tree text view with")
print(p.tree)
input()

print("Note!")
input()
print("The indexes in the tree represent the indexes in the sequences of `blast_seqs` but they are represented with indexes for clarity, to get names of the sequences for a given index, simply do: 'print(blast_seqs[0].name)'")
input()

# In[30]:


print(blast_seqs[0].name)
sep()

print("Graphs Section")
print("To visualize the graphs from the distance matrices, one only has to call the `graph(cut)` function for a given matrix, the `cut` parameter is a value grater than the maximum allowed distance for an edge to exist.")
input()

# In[31]:

print("Using the simple example for the matrix")
print(mxs)
input()


# In[32]:

print("Let's see the graph using a cut of 7")
input()
mxs.graph(cut=7)


# In[39]:

print("Now, Let's see the graph using a cut of 4 (less connections are expected)")
input()
mxs.graph(cut=4)


# In[33]:
sep()

print("For the Blast sequences we have the matrix")
print(mx)
input()


# In[38]:

print("Let's see the graph using a cut of 40")
mx.graph(40)

# In[34]:


print("Let's see the graph using a cut of 20")
mx.graph(20)


# In[35]:



print("Let's see the graph using a cut of 10")
mx.graph(10)



sep()

print("Now all the developed tests will be invoked, they thouroughly test the module and guarantee 100% coverage. After that we will also generate the html coverage, you can go to 'htmlcov/' and open index.html to see it. ok?")
input()
subprocess.call("coverage run -m unittest discover -v".split())
subprocess.call("coverage html".split())
print("Finally, the coverage report in the command line:")
subprocess.call("coverage report".split())
