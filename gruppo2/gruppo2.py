from Bio import Phylo, SeqIO
from io import StringIO
import subprocess as sp
import os
import itertools as it
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def distPlot(df):
    plt.hist(df[2], bins=100) 
    # questa funzione crea un istogramma delle distanze

protein_filepath = "cl3mente/PlantPath2023/data/data/protein" # il percorso del file in input con le sequenze FASTA
orthofinder = "OrthoFinder/orthofinder" # il percorso dove è installato orthofinder
# speciesNames = [x.lstrip(".") for x in os.listdir(protein_filepath)] # lista dei nomi dei file nella cartella protein

speciesNames = ["Mycoplasma_genitalium", "Mycoplasma_gallisepticum", "Mycoplasma_hyopneumoniae", "Mycoplasma_agalactiae"]
# prova Temporaneo 


"""
"for f in *fa ; do python ~/orthofinder_tutorial/OrthoFinder/tools/primary_transcript.py $f ; done" 
#comando preso dalla guida di orthofinder per velocizzare la run. devo usare primary_transcript.py che è nella cartella orthofinder
"""

#######################################################################################################################################
# ORTHOFINDER RUN                                                                                                                     #
# da aggiustare                                                                                                                       #
#######################################################################################################################################

"""
p1 = sp.Popen([orthofinder, "-f", protein_filepath]) # QUI inizia il comando che fa partire orthofinder
p1.wait() # aspetta che il comando finisca
if p1.returncode == 0:
    print("OrthoFinder run successful") # se il comando è andato a buon fine stampa questo messaggio
"""
    

########################################################################################################################################
# creating distMatrix                                                                                                                  #
# retrieving distances                                                                                                                 #
########################################################################################################################################

GeneTreesPath = "C:\Orthofinder\Results_Oct20\Gene_Trees" #directory Gene_Trees dall'output di Orthofinder
#TEMPORANEO: da aggiustare in modo che funzioni con il flusso di dati output in arrivo da Orthofinder

# getDistMatrix() restituisce un df con tutte le distanze tra i geni
# input: la directory in cui si trovano i file gene_tree.txt
def getDistMatrix(GeneTreesPath):

    # creo `trees`: contiene tutte le directory (path) dei file gene_tree.txt
    trees = []
    for path, names, files in os.walk(GeneTreesPath):
        for genetree in files:
            path_2 = os.path.join(path + "\\" + genetree)
            trees.append(path_2) 
        
    # creo `distances`: legge ciascun gene_tree.txt e calcola le distanze tra geni
    distances = []
    for i in trees:
        genTree = Phylo.read(i, "newick")
        clades = list(genTree.get_terminals())
        combs = list(it.combinations(clades, 2))
        
        for comb in combs:
            name1 = comb[0].name
            name2 = comb[1].name
            dist = genTree.distance(name1, name2)
            distances.append((name1, name2, dist))

    # `distances` è una lista di tuples: (stringa name1, stringa name2, float dist)

    # creo `distMatrix`: un dataframe che contiene tutte le distanze
    # in ogni riga di `distMatrix` c'è una coppia di geni omologhi e la loro distanza
    # creo un dataframe da `distances`
    distMatrix = pd.DataFrame(distances,
        columns = ["gene1", "gene2", "dist"])

    """ 
    distmatrix ha questo aspetto:

    gene1               gene2                   dist
    geneNameFromSp1     geneNameFromSp2         1.234 
    """

    return distMatrix


# splitMatrix() è una funzione che divide distMatrix in diversi df
# input: distMatrix (che è anche l'output di getDistMatrix())
# ciascun df contiene le distanze tra gli omologhi solo di due specie
# splitMatrix restituisce una lista di dataframe 
# (chiama "comparisons", perchè sono i paragoni di solo due specie)
# INFO: in questa funzione si usa `speciesNames`, da implementare
# serve una funzione tipo "getSpNames()" che legga i nomi dei file .fa
# e restituisce una lista di stringhe con i nomi delle specie
def splitMatrix(distMatrix):
    comparisons = []
    for comb in list(it.combinations(speciesNames, 2)): 
        #  per ogni combinazione di 2 specie 
        # (in nomi li prendo dai file stessi. da aggiustare 
        # (attualmente li prende perchè glieli ho scritti io in cima al codice))

        # crea "comparison" un subset del dataframe `distMatrix` che contiene 
        # solo le righe che iniziano con i nomi delle 2 specie
        # quindi contiene solo le distanze tra i geni delle 2 specie
        comparison = distMatrix[distMatrix["gene1"].str.startswith(comb[0]) & distMatrix["gene2"].str.startswith(comb[1])]
        
        # dagli un attributo ".name" che funziona così: "specie1 vs specie2"
        comparison.name = str(comb[0]) + " vs " + str(comb[1])

        # aggiungi `comparison` alla lista `comparisons`
        comparisons.append(comparison)
        # `comparisons` è una lista di dataframe, uno per ogni coppia di specie
        # ciascuna riga del df è una coppia di geni, con la distanza tra i 2 geni
    return comparisons


############################
#  finding HGT candidates  #
############################

# getHGT() è una funzione che prende in input un dataframe di distanze tra omologhi
# e calcola il threshold per i candidati HGT
# output: un df con i candidati HGT, nomi e distanze
def getHGT(comp):

    print(f"analyzing {comp.name}...")

    # grafico delle distanze
    # distPlot(comp)

    # calculating the HGT threshold and storing in min_dist
    # HGT threshold = media delle distanze - 3*deviazione standard
    # quindi sono i valori di distanza inferiori al confidence interval
    mean = np.mean(comp["dist"])
    std = np.std(comp["dist"])
    if mean < 3*std:
        min_dist = 0
        # qui potrei anche interrompere il loop: 
        # se il threshold è zero vuoldire che nessun gene è stato trasferito

        """ print(f"No HGT candidates found for {comp.name}")
        continue """
    else:
        min_dist = mean - 3*std

    # prendo tutti i geni sotto questa soglia
    # li metto in `candidates`
    candidates = comp[comp["dist"] < min_dist]

    # se trovi dei candidati, stampali
    if candidates.empty:
        print(f"No HGT candidates found for {comp.name}")
    else:
        return candidates


if __name__ == "__main__":
    distMatrix = getDistMatrix(GeneTreesPath)
    for c in splitMatrix(distMatrix):
        getHGT(c)