from Bio import Phylo, SeqIO
from io import StringIO
import subprocess as sp
import os
import itertools
import pandas as pd

protein_filepath = "data/data/protein" # il percorso del file in input con le sequenze FASTA
orthofinder = "OrthoFinder/orthofinder" # il percorso dove è installato orthofinder


"""
"for f in *fa ; do python ~/orthofinder_tutorial/OrthoFinder/tools/primary_transcript.py $f ; done" 
#comando preso dalla guida di orthofinder per velocizzare la run. devo usare primary_transcript.py che è nella cartella orthofinder
"""

############################
#     ORTHOFINDER RUN      #
############################

"""
p1 = sp.Popen([orthofinder, "-f", protein_filepath]) # QUI inizia il comando che fa partire orthofinder
p1.wait() # aspetta che il comando finisca
if p1.returncode == 0:
    print("OrthoFinder run successful") # se il comando è andato a buon fine stampa questo messaggio
"""
    

############################
#   retrieving distances   #
############################



GeneTreesPath = "C:\Orthofinder\Results_Oct20\Gene_Trees" #directory Gene_Trees dall'output di Orthofinder

trees = []
for path, names, files in os.walk(GeneTreesPath):
    for genetree in files:
        path_2 = os.path.join(path + "\\" + genetree)
        trees.append(path_2)

#in `trees` ci sono tutti i file Gene Tree
        

distances = []
for i in trees:
    genTree = Phylo.read(i, "newick")
    clades = list(genTree.get_terminals())
    combs = list(itertools.combinations(clades, 2))
    
    for comb in combs:
        name1 = comb[0].name
        name2 = comb[1].name
        dist = genTree.distance(name1, name2)
        distances.append((name1, name2, dist))



#in `distances` ci sono tutte le distanze tra le sequenze
distMatrix = pd.DataFrame(distances)
print(distMatrix)