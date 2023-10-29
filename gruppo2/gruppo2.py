from Bio import Phylo, SeqIO
from io import StringIO
import subprocess as sp
import os
import itertools as it
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import binascii

inputPath = "/mnt/c/Orthofinder/faino_prova_proteomi/protein" # il percorso del file in input con le sequenze FASTA
ofPath = "/OrthoFinder/orthofinder" # il percorso dove è installato orthofinder

def distPlot(df):
    plt.hist(df[2], bins=100) 
    # questa funzione crea un istogramma delle distanze

#######################################################################################################################################
# ORTHOFINDER RUN                                                                                                                     #
# da aggiustare                                                                                                                       #
#######################################################################################################################################

# generate a random code for a new folder in which orthofinder will store the results
def orthoResults(inputPath, ofPath):

    # generate a random hex code
    randomCode = binascii.hexlify(os.urandom(3)).decode('utf8')
    # assign it to the folder var, later pass it to the orthofinder command
    resultsFolder = inputPath + '/' + randomCode

    p1 = sp.Popen(["wsl ", ofPath, " -f ", inputPath, " -o ", resultsFolder], stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    output, errors = p1.communicate()

    p1.wait()  # Wait for the command to finish

    if errors:
        print(f"Orthofinder run failed. Error message: \n{errors}")
    if p1.returncode == 0:
        print("OrthoFinder run successful")
    return resultsFolder

########################################################################################################################################
# creating distMatrix                                                                                                                  #
# retrieving distances                                                                                                                 #
########################################################################################################################################

# getSpNames() è una funzione che legge il file SpeciesTree_rooted.txt
def getSpNames(ResultsPath):
    with open(ResultsPath + "\Species_Tree\SpeciesTree_rooted.txt", "r") as f:

        speciesTree = f.read()
        speciesTree = Phylo.read(StringIO(speciesTree), "newick")
        speciesNames = [clade.name for clade in speciesTree.get_terminals()]
        
        # speciesNames è una lista di stringhe con i nomi delle specie
        # speciesNames = ["Mycoplasma_genitalium", "Mycoplasma_gallisepticum", "Mycoplasma_hyopneumoniae", "Mycoplasma_agalactiae"]

        return speciesNames
    
# getDistMatrix() restituisce un df con tutte le distanze tra i geni
# input: la directory in cui si trovano i file gene_tree.txt
def getDistMatrix(Path):

    # creo `trees`: contiene tutte le directory (path) dei file gene_tree.txt
    trees = []
    for path, names, files in os.walk(Path):
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
def splitMatrix(distMatrix, ResultPath):
    comparisons = []
    for comb in list(it.combinations(getSpNames(ResultPath), 2)): 
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



    ResultsPath = orthoResults(inputPath, ofPath)

    """
    "for f in *fa ; do python ~/orthofinder_tutorial/OrthoFinder/tools/primary_transcript.py $f ; done" 
    # comando preso dalla guida di orthofinder per velocizzare la run. devo usare primary_transcript.py che è nella cartella orthofinder
    """
    #TEMPORANEO, da aggiustare con l'output della OF run
    GeneTreesPath = ResultsPath + "\Gene_Trees" 
    SpeciesTreesPath = ResultsPath + "\Species_Tree" #directory Species_Trees dall'output di Orthofinder

    distMatrix = getDistMatrix(ResultsPath + "\Gene_Trees")
    for c in splitMatrix(distMatrix, ResultsPath):
        getHGT(c)
