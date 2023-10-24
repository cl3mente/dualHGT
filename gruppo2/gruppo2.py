from Bio import Phylo, SeqIO
from io import StringIO
import subprocess as sp
import os

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

"""
for path, names, files in os.walk(orthofinder):
    if "SpeciesTree_rooted.txt" in files:
        tree = os.path.join(path, "SpeciesTree_rooted.txt") # il percorso dell'output di OrthoFinder con l'albero filogenetico

"""
if True: # un esempio di Newick phylogenetic tree per fare una prova
    
    #L'output di orthofinder è una stringa come questa dove sono codificati i nodi dell'albero e le bootstrap distances
    tree = "(Mycoplasma_hyopneumoniae:0.235223,(Mycoplasma_agalactiae:0.457015,(Mycoplasma_genitalium:0.453747,Mycoplasma_gallisepticum:0.414629)0.958955:0.381412)1:0.235223);"
    tree = StringIO(tree)


genTree = Phylo.read(tree, "newick")
Phylo.draw_ascii(genTree)

for clade in genTree.depths():
    print(clade.name, clade.branch_length)
