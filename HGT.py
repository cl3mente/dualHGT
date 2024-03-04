from pathlib import Path
from itertools import combinations
from collections import defaultdict
import plotly.express as px
import argparse
from Bio import Phylo, SeqIO
from io import StringIO
import subprocess as sp
import os
import itertools as it
import binascii
import os.path
import multiprocessing as mp
import tqdm
import pandas as pd
import subprocess
from matplotlib_venn import venn3
import numpy as np
import shutil
import matplotlib.pyplot as plt
import ete3

PARAAT = "ParaAT.pl  -h  %s  -a  %s  -n  %s   -p  %s  -o  %s -f axt"
KAKS = "KaKs_Calculator  -i %s -o %s -m %s"
ORTHOFINDER = "/data/bioinf2023/PlantPath2023/OrthoFinder/orthofinder -f %s -t %s -o %s %s"
# PARAAT = "ParaAT.pl  -h  %s  -a  %s  -n  %s  -m muscle -p  %s  -o  %s -f axt"


def arguments():
    parser = argparse.ArgumentParser(prog='HGT', description='This tool identify potential horizontally transfered genes between species')
    parser.add_argument('-i','--input',
                        help="In this directory should be present both reference and annotation with the same root", required= True)
    # parser.add_argument('references')
    parser.add_argument('-OF', '--orthofinder' ,type = str, default="")
    parser.add_argument('-OFr', '--orthofinderResults', type=str)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-nt', '--numberThreads',type=int, help='Number of parallel threads to use')
    parser.add_argument('-o', '--output', default = "output")
    args = parser.parse_args()
    return(args)

def gffread(path):
    """
    Reads GFF files, modifies them, and performs file operations.

    Parameters:
    -----------
        path (str): The path to the directory containing GFF and FASTA files from the species of interest.

    Returns:
    --------
        tuple: A tuple containing the paths to the protein files, the directory of two .fasta files (AA and DNA)
        and a gene association dictionary (matching a gene code to its original ID in the genome).
    """

    gene_association = {}   # gene_association is a dict used to save the mapping between the tag and the actual gene name

    res_path = os.path.join(path , "results")
    prot_path = os.path.join(res_path , "prot")
    cds_path =os.path.join(res_path , "cds")
    folder = Path(path)

    # from the specified directory (parameter 'path') retrieve the name of all .gff files
    gff = [".gff", ".gff3", ".gtf"]
    filename_gff = [os.path.join(folder, file.name) for file in folder.iterdir() if file.suffix in gff]
    count = 0
    species_name = []
    
    for file in filename_gff:
        with open(file + "_mod_gff", "w") as fho:
            with open(file, "r") as fh:
                for line in fh:
                    if len(line.split("\t")) == 9 and line.split("\t")[2].startswith("mRNA"):
                        count += 1
                        genes_collection = [file.rsplit(".", 1)[0].rsplit("/", 1)[1], line.replace("\t", " ").rstrip()]
                        line = line.rstrip() + ";HGT=gene" + str(count) #TODO to add an index for each genome file # ?
                        fho.write(line + "\n")
                        gene_association["gene" + str(count)] = genes_collection
                    else:
                        fho.write(line)

    extensions = [".fna", ".fasta", ".gff_mod_gff", ".gff3_mod_gff", ".gtf_mod_gff"]    # pick the modified gff files and the fasta files
    filename = [file.name for file in folder.iterdir() if file.suffix in extensions]
    species_unique = list(set(species_name))
    species_unique_number = [species_unique , len(species_unique)]

    if len(filename) % 2 != 0:
        print("Error: odd number of file, file do not match")
        raise SystemExit

    try:
        if not os.path.exists(res_path):
            Path(res_path).mkdir()
        if not os.path.exists(prot_path):
            Path(prot_path).mkdir()
        if not os.path.exists(cds_path):
            Path(cds_path).mkdir()
    except FileExistsError:
        print("Warning: res folder already exist")
    except Exception as e:
        print(f"Error creating directory '{res_path}': {e}")
        raise SystemExit

    # by filtering the data trough dictionary key we can run only 1 for cycle to create a nested list of file with the same name by using as a key the file name without extension

    file_dic = defaultdict(list)
    for file_name in filename:
        name, ext = file_name.rsplit(".", 1)
        file_dic[name].append(file_name)
    file_matched = [match for match in file_dic.values() if len(match) > 1]
    # we created a nested list with NOT sorted file that are matched

    for check in file_matched:
        if len(check) != 2:
            print("Error: file do not match")
            raise SystemExit

    for x in file_matched:

        # in these list comprehension we take the first value from the list that
        # that contains only 1 value that is the file name that ends with the specified
        # extensions then we save the file name by removing the extension

        gff_f_name = str(next((f_name for f_name in x if f_name.endswith(("_mod_gff"))), None))
        fna_f_name = str(next((f_name for f_name in x if f_name.endswith((".fna", ".fasta"))), None))
        res_name = str(gff_f_name.rsplit(".", 1)[0])

        # we specify where the files are going to be saved and the name of the saved file

        res_cds_file = cds_path + "/" + res_name + "_cds.fas"
        res_prot_file = prot_path + "/" + res_name + "_prot.faa"
        fasta_file = path + "/" + fna_f_name
        gff_file = path + "/" + gff_f_name

        # we prepared the command to be run by the command prompt as a list of strings
        # and then run the command, no output collection is needed as the command saves the files
        # in the specified directory

        argument_cds = ["gffread", "-w", res_cds_file, "-y", res_prot_file, "-F", "-g", fasta_file, gff_file]
        result_cds = subprocess.run(argument_cds)
        files = [res_cds_file, res_prot_file]


        for file in files:
            with open((file + "_mod.fasta"), "w") as Ffile:
                for record in SeqIO.parse(file, "fasta"):
                    tmp = record.description.split("HGT=")[1]
                    if ";" in tmp:
                        tmp = tmp.split(";")[0]
                    record.id = tmp
                    record.description = ""
                    SeqIO.write(record, Ffile, "fasta")

        os.remove(res_cds_file)
        os.remove(res_prot_file)

    prot_all_file =  res_path + '/proteinfilefinal.faa'
    cds_all_file =  res_path + '/cdsfilefinal.fas' 
    cmd = 'cat ' + prot_path + '/*.faa_mod.fasta' + ' > ' + prot_all_file   # this .fasta file is a collection of all the aminoacid sequences
    subprocess.run(cmd, shell=True)
    cmd = 'cat ' + cds_path + '/*.fas_mod.fasta' + ' > ' + cds_all_file # this .fasta file is a collection of all the coding DNA sequences
    subprocess.run(cmd, shell=True)

    # the output of the function is:
    #   the path to the two .fasta files (AA and DNA)
    #   the gene association dictionary
    #   the protein path

    return (prot_path, prot_all_file, cds_all_file, gene_association)

def orthoResults(inputPath, numberThreads, extra, verbose):
    """
    Runs Orthofinder on the specified directory and returns the path to the results folder.
    
    Args:
    -----
        `inputPath (str)`: The path to the directory containing the proteome files (.fasta format).
        `numberThreads (int)`: The number of threads to use for the Orthofinder run.
        `extra (str)`: Extra arguments to pass in the OrthoFinder command.
        `verbose (bool)`: If True, prints the output of the Orthofinder command.
        
    Returns:
    --------
        str: The path to the OrthoFinder results folder.
    """

    # generate a random hex code, to comply with multiple OrthoFinder runs in the same day and directory
    randomCode = binascii.hexlify(os.urandom(3)).decode('utf8')

    # assign it to the folder var, later pass it to the orthofinder command
    resultsFolder = inputPath + '/' + randomCode

    # write the linux command to run orthofinder
    runortho = ORTHOFINDER % (inputPath,str(numberThreads), resultsFolder, extra)

    # run the command with as a subprocess
    p1 = sp.Popen(runortho,shell=True)
    stdout, stderr = p1.communicate()
    if verbose:
        print(stdout)
        print(stderr)

    if p1.returncode != 0:
        print(f"OrthoFinder failed: error code {p1.returncode}")

    # retrieve the most recent folder created in the result directory, which will contain the OrthoFinder results, and return its path
    arr = os.listdir(resultsFolder)
    resultsFolder += "/" + arr[0]
    return resultsFolder

def getSpNames(ResultsPath):
    """
    Retrieve the names of the species from the species tree file, created by OrthoFinder.
    ## Args:
        ResultsPath (str): The path to the Orthofinder results folder.
    ## Returns:
        list: A list of strings containing the names of the species.
    """
    with open(ResultsPath + "/Species_Tree/SpeciesTree_rooted.txt", "r") as f:

        speciesTree = f.read()
        speciesTree = Phylo.read(StringIO(speciesTree), "newick")
        speciesNames = [clade.name for clade in speciesTree.get_terminals()]

        # speciesNames è una lista di stringhe con i nomi delle specie
        return speciesNames

def read_tree(info):
    """
    Reads the a genetree file and returns a list containing gene pairs, their orthogroup, and their distance.
    
    ## Args:
        info (tuple): A tuple containing the path to the genetree file and the species list.
    ## Returns:
        list: A list containing gene pairs and their distances. This list corresponds to a row of the final dataframe used to estimate HGT occurrence.
    
    This function is called in `parseOrthofinder()` and executed in parallel with mp.Pool()
    """

    genetree, species_list = info
    distances=[]
    genTree = Phylo.read(genetree, "newick")
    clades = list(genTree.get_terminals())
    combs = list(it.combinations(clades, 2))

    for comb in combs:

        name1 = comb[0].name
        name2 = comb[1].name

        if name1 < name2:
            name1, name2 = name2, name1

        dist = genTree.distance(name1, name2)
        #
        # sp1 = sorted([os.path.commonprefix([i, name1]) for i in species_list], key=len)
        # sp2 = sorted([os.path.commonprefix([i, name2]) for i in species_list], key=len)
        #
        # species_comb = sp1[-1] + "_vs_" + sp2[-1]

        type = "tree"
        name1 = name1.split("_")[-1]
        name2 = name2.split("_")[-1]
        distances.append([name1, name2, float(dist), genetree.split("/")[-1].split("_")[0], type])
        # the final output is `distances`, a list that looks like this:
        # [["gene1", "gene2", "dist", "OG", "type"]]

    return distances

def parseOrthofinder(ResultsPath: str, threads: int):

    """
    A function to parse the Orthofinder results and return a list of entries containing the gene pairs and their distances.
    
    ## Args:
        ResultsPath (str): The path to the Orthofinder results folder.
        threads (int): The number of threads to use for the parsing.
    ## Returns:
        list: Data containing the gene pairs and their distances. The data comes in the form of a list of lists.
        each list is to be considered as a data entry and will contain these columns:
            `gene1` | `gene2` | `dist` inferred by OrthoFinder | `OG` | `type` = "tree" 
    """
    
    # retrieve the species names from the species tree file
    species_list = getSpNames(ResultsPath)

    # retrieve the gene trees from the gene trees folder
    GeneTreesPath = ResultsPath + "/Gene_Trees"
    for (dir_path, dir_names, file_names) in os.walk(GeneTreesPath):
        files = file_names
    tree_abs_path = [(os.path.join(GeneTreesPath, file), species_list) for file in files]
    
    distances = []
    # read the gene trees in parallel with multiprocessing.Pool()
    with mp.Pool(threads) as p:
        for x in tqdm.tqdm(p.imap_unordered(read_tree, tree_abs_path), total=len(tree_abs_path), desc="Reading GeneTrees..."):
            distances.append(x) # appends each entry to `distances`
            pass
    
    distances = [item for sublist in distances for item in sublist] # turn it into a flat list
    
    return distances

def kaksparallel(file):
    """
    Runs the KaKs calculator program in the shell.
    
    ## Args:
        file (str): The path to the .axt file to be processed.
    ## Returns:
        list: A list containing the gene pair and their distance.

    References:
    -----------
    Nei, M., & Gojobori, T. (1986). Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. Molecular biology and evolution, 3(5), 418–426. https://doi.org/10.1093/oxfordjournals.molbev.a040410
    """

    # run the 'KAKS' command in the shell
    output = file + ".kaks"
    if not os.path.exists(output) or not os.path.getsize(output) > 0:
        runkaks = KAKS % (file, #the .axt file passed as input
                          output,   # the output file
                          "NG") # NG is the model used for the calculation
        run = subprocess.Popen(runkaks, shell=True,stderr=subprocess.PIPE ,stdout=subprocess.PIPE)
        out, err = run.communicate()
        # print(out)
        print(err)

    # read the output file and return the gene pair and their distance
    my_list = read_kaks_file(output)

    return my_list

def read_kaks_file(kaks_filename):
    with open(kaks_filename, newline='') as resultskaks:
        next(resultskaks)
        for line in resultskaks:
            if line.strip():
                # take the first and fourth elements of the file
                # corresponding to the sequence name and the Ka/Ks ratio
                my_list = [line.split('\t')[i] for i in [0, 3]]
    return my_list

def parseKaKs(ResultsPath, threads, proteinfilefinal,cdsfilefinal):
    """
    A function to combine the gene pairs from the Orthofinder results and the KaKs distances.
    ## Args:
        ResultsPath (str): The path to the Orthofinder results folder.
        threads (int): The number of threads to use for the parsing.
        proteinfilefinal (str): The path to the protein file.
        cdsfilefinal (str): The path to the CDS file.
    ## Returns:
        list: A list containing the gene pairs and their distances. The list will contain this data:
            `gene1` | `gene2` | `dist` = Ka/Ks ratio | `OG` | `type` = "kaks"
    """
    fileThreads = ResultsPath + "/proc"
    with open(fileThreads, "w") as fh:
        fh.write(str(threads) +"\n")

    # create a dictionary to match genes to orthogroups, output
    data = pd.read_csv(ResultsPath  + "/Orthogroups/Orthogroups.tsv", sep="\t")
    data_dict = {}
    data.fillna('empty', inplace = True)

    for row in data.itertuples(index = False):
        values=[]
        for i in range(1, len(data.columns)):
            key = row[0]
            values.append(row[i].split(', '))
            data_dict[key] = values

    # a dictionary that contains gene pairs matched to their orthogroup.
    # questo va cambiato per evitare la run di paraAT
    dict_match = {}
    file_out = os.path.join('/tmp/output.txt')
    with open(file_out, 'w') as output_file:
        for group_name, group_values in data_dict.items():
            for genes1, genes2 in combinations(group_values, 2):
                for gene1 in genes1:
                    if gene1 == "empty":
                        break
                    for gene2 in genes2:
                        if gene2 == "empty":
                            break
                        #if gene1 > gene2:
                        #    gene1, gene2 = gene2, gene1 # impose an order to avoid duplication of pairs
                        line = f"{gene1}\t{gene2}\n"
                        output_file.write(line)
                            # data comes in this format:
                            #   {"gene1-gene2": "OG"}
                        dict_match[gene1 + "-" + gene2] = group_name

    kaksfolder = ResultsPath + "/kaksfolder"

    # check if there's any .axt filepath from paraAT

    axtFiles = []
    for filename in os.listdir(kaksfolder):
        if filename.endswith('.axt'):
            axtFiles.append(os.path.join(kaksfolder, filename))


    if axtFiles:
        pass
    else:
        # run paraAT in the shell.
        # This will create .axt files from the protein and CDS files
        runparaAT = PARAAT % (file_out, proteinfilefinal, cdsfilefinal, fileThreads, kaksfolder)
        print(runparaAT)
        run = subprocess.Popen(runparaAT, shell=True, cwd=ResultsPath,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = run.communicate()
        if arg.verbose:
            print(out)
            print(err)
        for filename in os.listdir(kaksfolder):
            if filename.endswith('.axt'):
                axtFiles.append(os.path.join(kaksfolder, filename))

    var = []

    kaks_filepaths = []
    for filename in os.listdir(kaksfolder):
        if filename.endswith('.kaks'):
            kaks_filepaths.append(os.path.join(kaksfolder, filename))

    if kaks_filepaths:
        for file in kaks_filepaths:
            row = read_kaks_file(file)
            var.append(row)
    else: # run 'kaksparallel' function on the .axt files in parallel with multiprocessing.Pool()
        with mp.Pool(threads) as p:
            with tqdm.tqdm(total=len(axtFiles), desc="running kaks...") as pbar:
                for x in p.imap_unordered(kaksparallel, axtFiles):
                    # x is a list that looks like this:
                    #   ['seq_(pair?)_name', 'Ks']
                    var.append(x)
                    # therefore var will be a list of lists:
                    #   [
                    #       ['seq_(pair?)_name1', 'Ks'],
                    #       ['seq_(pair?)_name2', 'Ks'],
                    #       ...
                    #   ]
                    pbar.update()

    # create an 'info' list with KaKs score and the gene pairs.
    #   data comes in this format:
    #   [[gene1, gene2, OG, Ks]]
                
    distances = []

    for entry in var:
        OG = dict_match[entry[0]] # retrieve orthogroup by matching with the seq_pair name
        genes = entry[0].rsplit("-") # retrieve the original gene names

        del entry[0] # delete the seq_pair name

        row = genes + entry + [OG, "kaks"]
        distances.append(row)

    return distances

def append_species(entry_list, dict_species):
    """
    Transform an entry list (a list of rows) into a canonical pd.DataFrame.
    Add to the dataframe the information regarding the species to which the genes examined belong to
    
    Args:
    -----
    entry_list: list
        A list of lists containing the gene pairs and their distances.
    dict_species: dict
        A dictionary containing the gene names and the species to which they belong.
    
    Returns:
    --------
    pd.DataFrame
        A dataframe containing the gene pairs and their distances, with the species information added.
    """

    final = []
    for line in entry_list:
        species = []
        if dict_species[line[0]]:
            species.append(dict_species[line[0]][0])
        if dict_species[line[1]]:
            species.append(dict_species[line[1]][0])
        species.sort()
        line.append("_".join(species))
        final.append(line)


    matrix = pd.DataFrame(final,
         columns = ["gene_1", "gene_2", "dist", "OG",  "type", "species"])

    return matrix

def getMeanDist(comp):

    print(f"analyzing {comp.iloc[0,4]}...")

    mean = comp['dist'].mean()
    comp["mean_dist"] = mean - comp["dist"]
    comp.sort_values(by=["mean_dist"], inplace=True, ascending=False)

    return comp[comp.mean_dist > 0]

def getHGT(matrix, dict_species):
    """
    Main function to identify potential HGT events between species.
    
    Parameters:
    -----------
    `matrix` : pd.DataFrame
        A list of dataframes containing the gene pairs and their distances. The dataframes will be derived from the outputs of KaKs and OrthoFinder.

    Returns:
    --------
        `matrix` : pd.DataFrame
            A dataframe with the `HGT` column added, containing the HGT score for each gene pair.
    """

    matrix2 = append_species(matrix, dict_species)

    # format the pd.DataFrame from list
    matrix2 = matrix2[matrix2['dist'] != "NA"]  # remove NAs
    matrix2['dist'] = pd.to_numeric(matrix2['dist'], downcast='float')
    matrix2['HGT'] = None   # initialize an empty column for the HGT score
    sp_pairs = matrix2['species'].unique()    

    for sp in sp_pairs:
        threshold = matrix2[matrix2['species'] == sp]['dist'].quantile(.05) # get the 5th percentile of the distance.
        matrix2.loc[matrix2['species'] == sp, 'HGT'] = matrix2['dist'] < threshold # set the 'HGT' variable to True if the distance is less than the threshold
    
    return matrix2

    """ 
    for key in new_dfs: 
        # for each species pair create a dataframe
        dataFrame = new_dfs[key]
        threshold = dataFrame['dist'].quantile(.05) # get the 5th percentile of the distance.

        # in the species dataframe, create the 'set' column 
        # the 'set' variable is informative of whether a candidate HGT event
        dataFrame["HGT"] = dataFrame['dist' < threshold]
        distMatrix3 = dataFrame[dataFrame['dist'] < threshold]

        list_value =distMatrix3.values.tolist()

        if key in dict_dataframes:
            dict_dataframes[key] = dict_dataframes[key] + [list_value]
        else:
            dict_dataframes[key] = [list_value]

        data = dataFrame[["OG", "HGT"]]
        data = data.sort_values(by=["HGT"])
        data["HGT"] = data["HGT"].astype(int)
        dataconcatenate.concat(dataFrame)

        true_false_dict = data.set_index("OG")["HGT"].to_dict() 

    kaks_table= []
    tree_table = []
    for key in dict_dataframes:
        if len(dict_dataframes[key]) == 2:
            species1 = dict_dataframes[key][0]
            species2 = dict_dataframes[key][1]
            for line in species1:
                for res in species2:
                    if line[0] in res:
                        kaks_table.append(line)
                        tree_table.append(res)
    """
    """
     dist_ks = pd.DataFrame(kaks_table,
                               columns=["gene_1", "gene_2", "dist", "OG", "type", "species"])
    dist_ks_sorted = dist_ks.sort_values(by='dist')

    dist_tree = pd.DataFrame(tree_table,
                               columns=["gene_1", "gene_2", "dist", "OG", "type", "species"])
    dist_tree_sorted = dist_tree.sort_values(by='dist') 
    """

def get_topology(ResultsPath):
    """
    Parameters
    ----------
    ResultsPath : str
        The path to the Orthofinder results folder.

    Returns
    ---------
    list:
        A list of orthogroups with significantly different topology from that of the average species tree.
    """
    # get the species tree as a reference
    single_tree_folder = os.path.join(ResultsPath, "Gene_Trees/")
    species_tree = os.path.join(ResultsPath, "Species_Tree/SpeciesTree_rooted.txt")

    with open(species_tree, "r") as fh:
        spRoot_single = fh.read()
    species_tree_ete = ete3.Tree(spRoot_single)
    for node in species_tree_ete.traverse():
        if node.is_leaf():
            node.name = node.name.replace(".", "_")
    print(species_tree_ete)

    dict_topology = {} # a dictionary with OG as key and their 'HGT score' (1 or 0) as value
    list_keep = []

    for root, dir, files in os.walk(single_tree_folder):
        for file in files:
            og = file.split("_")[0]
            dict_topology[og] = 0
            file_og = os.path.join(single_tree_folder, file)
            with open(file_og, "r") as fh :
                og_single = fh.read()
                og_tree = ete3.Tree(og_single)
                for node in og_tree.traverse():
                    if node.is_leaf():
                        node.name = node.name.split("_gene")[0]
                try:
                    diff = og_tree.compare(species_tree_ete)    # use the 'compare()' function from ete3 to compute the Robinson-Foulds distance
                    if diff["rf"] > 0:  # if the distance from the average species tree is greater than 0, store the 
                        list_keep.append(og)
                        dict_topology[og] = 1
                except Exception as x:
                    list_keep.append(og)
                    dict_topology[og] = 1
                    continue

    list_uniq = list(set(list_keep))    # remove duplicates from the list of HGT candidate OGs

    return list_uniq

def vennPlot(kaks_OG_list, tree_OG_list, topology_OG_list):
    """

    Make a Venn Diagram. Easy.

    Parameters:
    ----------
    dist_matrix_kaks : pd.DataFrame
        A list of with the Ks distances sorted in ascending order.
    dist_matrix_tree : pd.DataFrame
        A dataframe with the distance inferred by OrthoFinder sorted in ascending order.
    list_topology : list
        A list of orthogroups with significantly different topology from that of the average species tree.

    """
    plt.figure(figsize=(4, 4))
    venn3([set(kaks_OG_list),
           set(tree_OG_list),
           set(topology_OG_list)],
          set_labels=('KaKs', 'Tree', 'Topology')
          )
    plt.show()

def plotData(df):
    """
    A final function to plot the data.

    Parameters:
    ----------
    final_dataset : pd.DataFrame
        A dataframe with all information to plot:
            gene names, Ka/Ks ratio and OrthoFinder distance, species identity
    """

    fig = px.violin(df,
                    y = "dist", 
                    x = "species", 
                    color = "type", 
                    box = True, 
                    hover_data = df.columns)
    fig.show()

if __name__ == "__main__":

    arg = arguments()

    prot_path,prot_all_file, cds_all_file,dict_species = gffread(arg.input)

    if not arg.orthofinderResults:
        ResultsPath = orthoResults(prot_path, arg.numberThreads,arg.orthofinder,arg.verbose)
    else:
        ResultsPath = arg.orthofinderResults

    # shutil.rmtree("/data/bioinf2023/PlantPath2023/data/genomeANDgff/results/prot/16cb6a/Results_Nov28/kaksfolder", )
    # ResultsPath = "/data/bioinf2023/PlantPath2023/data/genomeANDgff/results/prot/16cb6a/Results_Nov28"
    # prot_all_file = "/data/bioinf2023/PlantPath2023/data/genomeANDgff/results/proteinfilefinal.faa"
    # cds_all_file = "/data/bioinf2023/PlantPath2023/data/genomeANDgff/results/cdsfilefinal.fas"

    dist_matrix_tree = parseOrthofinder(ResultsPath, arg.numberThreads)
    dist_matrix_tree = getHGT(dist_matrix_tree, dict_species)

    dist_matrix_kaks = parseKaKs(ResultsPath,arg.numberThreads,prot_all_file,cds_all_file)
    dist_matrix_kaks = getHGT(dist_matrix_kaks, dict_species)

    list_kaks = dist_matrix_kaks.loc[dist_matrix_kaks['HGT'] == True,'OG'].to_list()
    list_tree = dist_matrix_tree.loc[dist_matrix_tree['HGT'] == True,'OG'].to_list()
    list_topology = get_topology(ResultsPath)

    # TODO modificare arg.verbose che va qua
    if False:
    #if arg.verbose: # show the topology of the candidate HGT orthogroups
        print("printing the topology of candidate HGT orthogroups:\n")
        for OG in list_topology:
            single_tree_folder = os.path.join(ResultsPath, "Gene_Trees/")
            file_og = os.path.join(single_tree_folder, OG + "_tree.txt")
            with open(file_og, "r") as fh :
                og_single = fh.read()
                og_tree = ete3.Tree(og_single)
                print(OG)
                print(og_tree)

    vennPlot(list_kaks, 
             list_tree, 
             list_topology)

    intersection = list(set([i for i in list_kaks if i in list_tree and i in list_topology]))
    print(len(intersection))
    inters_path = os.path.join("intersection.tsv")
    if os.path.exists(inters_path):
        os.remove(inters_path)
    with open(inters_path, 'x') as f:
        for i in intersection:
            f.write(i + '\n')

    single_tree_folder = os.path.join(ResultsPath, "Gene_Trees/")
    for i in intersection:
        file_og = os.path.join(single_tree_folder, i + "_tree.txt")
        with open(file_og, "r") as fh:
            og_single = fh.read()
        og_tree = ete3.Tree(og_single)
        print(i)
        print(og_tree)



    plot_matrix = pd.concat([dist_matrix_kaks, dist_matrix_tree], axis=0)

    full_matrix = pd.pivot(plot_matrix, values=['dist', 'HGT'], columns='type', index=['gene_1','gene_2', 'OG']).reset_index()
    full_matrix.columns = ['gene_1', 'gene_2', 'OG', 'dist_kaks', 'dist_tree', 'HGT_kaks', 'HGT_tree']

    full_matrix['HGT_topology'] = full_matrix['OG'].isin(list_topology).astype(int) # add the topology score if the OG appears in `list_topology`

    # compute the final HGT score and sort the dataframe
    full_matrix['HGT'] = full_matrix['HGT_kaks'] + full_matrix['HGT_tree'] + full_matrix['HGT_topology']
    full_matrix.sort_values(by=['HGT'], inplace=True, ascending=False)

    # TODO go back to the original protein name with dict_match

    match = {}
    for key, value in dict_species.items():
        mkey = key
        sp_name, genID = value[0], value[1]
        genID = genID.split('ID=')[-1]
        mvalue = '_'.join([sp_name,genID])
        match[mkey] = mvalue

    print('a')
    print(match['gene1'])

    full_matrix.replace({'gene_1':match, 'gene_2':match}, inplace=True)

    names = {}
    for i in full_matrix['gene_1']:

        if i in match:

            names.append(match[i])
    full_matrix['gene_1'] = names


    print('b')

    # write the final output to a .tsv file; in the form of a dataframe with scores from the KaKs, tree, and topology analyses
    with open('HGT_candidates.tsv', 'w') as f:
        full_matrix.to_csv(f, sep='\t', index=False)
    
    plotData(plot_matrix)