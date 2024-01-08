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
from matplotlib_venn import venn2
import shutil
import matplotlib.pyplot as plt
import ete3




PARAAT = "ParaAT.pl  -h  %s  -a  %s  -n  %s   -p  %s  -o  %s -f axt"
KAKS = "KaKs_Calculator  -i %s -o %s -m %s"
ORTHOFINDER = "/home/lfaino/Downloads/OrthoFinder2.3.8/orthofinder -f %s -t %s -o %s %s"
PARAAT = "ParaAT.pl  -h  %s  -a  %s  -n  %s  -m muscle -p  %s  -o  %s -f axt"


def arguments():
    parser = argparse.ArgumentParser(prog='HGT', description='This tool identify potential horizontally transfered genes between species')
    parser.add_argument('-i','--input',
                        help="In this directory should be present both reference and annotation with the same root", required= True)
    # parser.add_argument('references')
    parser.add_argument('-OF', '--orthofinder' ,type = str, default="")      # option that takes a value
    parser.add_argument('-OFr', '--orthofinderResults', type=str)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-nt', '--numberThreads',type=int, help='Number of parallel threads to use')
    parser.add_argument('-o', '--output', default = "output")
    args = parser.parse_args()
    return(args)

def kaksparallel(file):

    output = file + ".kaks"
    if not os.path.exists(output) or not os.path.getsize(output) > 0:
        runkaks = KAKS % (file, output,"NG")
        run = subprocess.Popen(runkaks, shell=True,stderr=subprocess.PIPE ,stdout=subprocess.PIPE)
        out, err = run.communicate()
        #print(out)
        print(err)

    with open(output, newline='') as resultskaks:
        next(resultskaks)
        for line in resultskaks:
            if line.strip():
                my_list = [line.split('\t')[i] for i in [0, 3]]

    return my_list

def geneComb (ResultsPath, threads, proteinfilefinal,cdsfilefinal):


    fileThreads = ResultsPath + "/proc"
    with open(fileThreads, "w") as fh:
        fh.write(str(threads) +"\n")

    data = pd.read_csv(ResultsPath  + "/Orthogroups/Orthogroups.tsv", sep="\t")
    data_dict = {}
    data.fillna('empty', inplace = True)

    for row in data.itertuples(index = False):
        values=[]
        for i in range(1, len(data.columns)):
            key = row[0]
            values.append(row[i].split(', '))
            data_dict[key] = values
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
                        line = f"{gene1}\t{gene2}\n"
                        output_file.write(line)
                        dict_match[gene1 + "-" + gene2] = group_name

    kaksfolder = ResultsPath + "/kaksfolder"
    # runpataAT = PARAAT % (file_out, proteinfilefinal, cdsfilefinal, "./proc", kaksfolder)
    # print(runpataAT)
    #
    # run = subprocess.Popen(runpataAT, shell=True, cwd=ResultsPath,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # out,err = run.communicate()

    axtFiles = []
    for filename in os.listdir(kaksfolder):
        if filename.endswith('.axt'):
            axtFiles.append(os.path.join(kaksfolder, filename))

    var = []
    with mp.Pool(threads) as p:
        with tqdm.tqdm(total=len(axtFiles), desc="running kaks...") as pbar:
            for x in p.imap_unordered(kaksparallel, axtFiles):
                var.append(x)
                pbar.update()

    # create dataframe to pass to group4 with results from kaks
    info = []
    for line in var:
        og = dict_match[line[0]]
        genes = line[0].rsplit("-")
        del line[0]
        line.append(og)
        supp = genes + line + ["kaks"]
        info.append(supp)

    # step3_results = pd.DataFrame(var, columns=["gene1", "gene2","dist", "OG","type"])

    return(info)

def gffread(path):
    gene_association = {}
    # from the specified directory we save the name of all file that ends with
    # a fasta extension or a gff extension
    # path = "/data/bioinf2023/PlantPath2023/data/genomeANDgff"
    res_path = os.path.join(path , "results")
    prot_path = os.path.join(res_path , "prot")
    cds_path =os.path.join(res_path , "cds")
    folder = Path(path)

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
                        line = line.rstrip() + ";HGT=gene" + str(count) #TODO to add an index for each genome file
                        fho.write(line + "\n")
                        gene_association["gene" + str(count)] = genes_collection
                    else:
                        fho.write(line)
    extensions = [".fna", ".fasta", ".gff_mod_gff", ".gff3_mod_gff", ".gtf_mod_gff"]
    filename = [file.name for file in folder.iterdir() if file.suffix in extensions]
    species_unique = list(set(species_name))
    species_unique_number = [species_unique , len(species_unique)]
    if len(filename) % 2 != 0:
        print("Error: odd number of file, file do not mach")
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

    # by filtering the data trough dictionary key we can run only 1 for cycle to
    # create a nested list of file with te same name by using as a key the
    # file name without extension
    file_dic = defaultdict(list)
    for file_name in filename:
        name, ext = file_name.rsplit(".", 1)
        file_dic[name].append(file_name)
    file_matched = [mach for mach in file_dic.values() if len(mach) > 1]
    # we created a nested list with NOT sorted file that are matched

    for ceck in file_matched:
        if len(ceck) != 2:
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
    cmd = 'cat ' + prot_path + '/*.faa_mod.fasta' + ' > ' + prot_all_file
    subprocess.run(cmd, shell=True)
    cmd = 'cat ' + cds_path + '/*.fas_mod.fasta' + ' > ' + cds_all_file
    subprocess.run(cmd, shell=True)

    return(prot_path,prot_all_file, cds_all_file, gene_association)

def plotData(final_dataset):

    final_df = pd.DataFrame(final_dataset)

    fig = px.violin(final_df,y = "dist", x = "species", color = "type", box = True, hover_data = final_df.columns)
    fig.show()

#######################################################################################################################################
# ORTHOFINDER RUN                                                                                                                     #
# da aggiustare                                                                                                                       #
#######################################################################################################################################

# generate a random code for a new folder in which orthofinder will store the results
def orthoResults(inputPath, numberThreads, extra, verbose):

    # generate a random hex code
    randomCode = binascii.hexlify(os.urandom(3)).decode('utf8')
    # assign it to the folder var, later pass it to the orthofinder command
    resultsFolder = inputPath + '/' + randomCode
    runortho = ORTHOFINDER % (inputPath,str(numberThreads), resultsFolder, extra)
    print(runortho)

    p1 = sp.Popen(runortho,shell=True)
    stdout, stderr = p1.communicate()
    if verbose:
        print(stdout)
        print(stderr)

    if p1.returncode != 0:
        print(f"OrthoFinder failed: error code {p1.returncode}")
    arr = os.listdir(resultsFolder)
    resultsFolder += "/" + arr[0]
    return resultsFolder

########################################################################################################################################
# creating distMatrix                                                                                                                  #
# retrieving distances                                                                                                                 #
########################################################################################################################################

# getSpNames() è una funzione che legge il file SpeciesTree_rooted.txt
def getSpNames(ResultsPath):
    with open(ResultsPath + "/Species_Tree/SpeciesTree_rooted.txt", "r") as f:

        speciesTree = f.read()
        speciesTree = Phylo.read(StringIO(speciesTree), "newick")
        speciesNames = [clade.name for clade in speciesTree.get_terminals()]

        # speciesNames è una lista di stringhe con i nomi delle specie
        return speciesNames

# read_tree() legge un genetree e dà una lista con  input: la directory in cui si trovano i file gene_tree.txt
def read_tree(info):

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

    return distances

def append_species(kaks_tree, dict_species):

    final = []
    for line in kaks_tree:
        species = []
        if dict_species[line[0]]:
            species.append(dict_species[line[0]][0])
        if dict_species[line[1]]:
            species.append(dict_species[line[1]][0])
        species.sort()
        line.append("_".join(species))
        final.append(line)


    matrix = pd.DataFrame(final,
         columns = ["gene1", "gene2", "dist", "OG",  "type", "species"])

    return matrix

def splitMatrix(distMatrix):

    comparisons = []
    unique_vals = distMatrix["species"].unique()
    print(len(unique_vals))
    for i in unique_vals:
        comb = pd.DataFrame(distMatrix[distMatrix["species"] == i])
        comparisons.append(comb)
    return comparisons

def getMeanDist(comp):

    print(f"analyzing {comp.iloc[0,4]}...")

    mean = comp['dist'].mean()
    comp["mean_dist"] = mean - comp["dist"]
    comp.sort_values(by=["mean_dist"], inplace=True, ascending=False)

    return comp[comp.mean_dist > 0]

def getHGT(Matrixs):

    dict_dataframes = {}
    og_TF = {"kaks":{},"tree":{}}
    test = {}
    for distMatrix in Matrixs:
        distMatrix1 = pd.DataFrame(distMatrix,
                                   columns = ["gene1", "gene2", "dist", "OG",  "type", "species"])
        distMatrix2 = distMatrix1[distMatrix1['dist'] != "NA"]
        distMatrix2['dist'] = pd.to_numeric(distMatrix2['dist'], downcast='float')
        groups = distMatrix2['species'].unique()
        new_dfs = {}
        for group in groups:
            name = str(group)
            new_dfs[name] = distMatrix2[distMatrix2['species'] == group]

        dataconcatenate = pd.DataFrame()
        for key in new_dfs:
            dataFrame = new_dfs[key]
            quantile = dataFrame['dist'].quantile(.05)
            dataFrame["set"] = dataFrame['dist'] < quantile
            distMatrix3 = dataFrame[dataFrame['dist'] < quantile]
            list_value =distMatrix3.values.tolist()
            if key in dict_dataframes:
                dict_dataframes[key] = dict_dataframes[key] + [list_value]
            else:
                dict_dataframes[key] = [list_value]
            data = dataFrame[["OG", "set"]]
            data = data.sort_values(by=["set"])
            data["set"] = data["set"].astype(int)
            dataconcatenate.concat(dataFrame)

            true_false_dict = data.set_index("OG")["set"].to_dict()

            # for key in true_false_dict["OG"]:
            #     if
            #     test[true_false_dict["OG"][key]] = true_false_dict["set"][key]
            #
            # for ind in dataFrame.index:
            #     if dataFrame['set'][ind]:
            #         print(dataFrame['OG'][ind])

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

    dist_ks = pd.DataFrame(kaks_table,
                               columns=["gene1", "gene2", "dist", "OG", "type", "species"])
    dist_ks_sorted = dist_ks.sort_values(by='dist')
    dist_tree = pd.DataFrame(tree_table,
                               columns=["gene1", "gene2", "dist", "OG", "type", "species"])
    dist_tree_sorted = dist_tree.sort_values(by='dist')
    return(dist_ks_sorted,dist_tree_sorted)

def venDiagram(comps):

    dist_matrix_kaks = comps[0]
    dist_matrix_tree = comps[1]



    # kaks = [line[3] for line in dist_matrix_kaks]
    # tree = [line[3] for line in dist_matrix_tree]

    plt.figure(figsize=(4, 4))

    venn2([set(dist_matrix_kaks['OG'].to_list()),
           set(dist_matrix_tree['OG'].to_list())],
          set_labels=('KaKs', 'Tree')
          )


    plt.show()
    return ()

def parseOrthofinder(ResultsPath, threads):
    species_list = getSpNames(ResultsPath)

    GeneTreesPath = ResultsPath + "/Gene_Trees"
    for (dir_path, dir_names, file_names) in os.walk(GeneTreesPath):
        files = file_names
    tree_abs_path = [(os.path.join(GeneTreesPath, file), species_list) for file in files]
    distances = []
    with mp.Pool(threads) as p:
        for x in tqdm.tqdm(p.imap_unordered(read_tree, tree_abs_path), total=len(tree_abs_path)):
            distances.append(x)
            pass
    distances = [item for sublist in distances for item in sublist]

    #matrix = dist_matrix(distances)

    return(distances)

def topology(ResultsPath):
    single_tree_folder = os.path.join(ResultsPath, "Gene_Trees/")
    species_tree = os.path.join(ResultsPath, "Species_Tree/SpeciesTree_rooted.txt")
    with open(species_tree, "r") as fh:
        spRoot_single = fh.read()
    species_tree_ete = ete3.Tree(spRoot_single)
    for node in species_tree_ete.traverse():
        if node.is_leaf():
            node.name = node.name.replace(".", "_")
    print(species_tree_ete)
    dict_topology = {}
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
                    diff = og_tree.compare(species_tree_ete)
                    if diff["rf"] > 0:
                        list_keep.append(og)
                        dict_topology[og] = 1
                except Exception as x:
                    list_keep.append(og)
                    dict_topology[og] = 1
                    continue
    list_uniq = list(set(list_keep))
    for tree in list_uniq:
        file_og = os.path.join(single_tree_folder, tree + "_tree.txt")
        with open(file_og, "r") as fh :
            og_single = fh.read()
            og_tree = ete3.Tree(og_single)
            print(tree)
            print(og_tree)

    return (list_uniq, dict_topology)

if __name__ == "__main__":
    arg = arguments()
    prot_path,prot_all_file, cds_all_file,dict_species = gffread(arg.input)
    if not arg.orthofinderResults:
        ResultsPath = orthoResults(prot_path, arg.numberThreads,arg.orthofinder,arg.verbose)
    else:
        ResultsPath = arg.orthofinderResults
    dist_matrix_tree =parseOrthofinder(ResultsPath, arg.numberThreads)
    # shutil.rmtree("/data/bioinf2023/PlantPath2023/data/genomeANDgff/results/prot/16cb6a/Results_Nov28/kaksfolder", )
    # ResultsPath = "/data/bioinf2023/PlantPath2023/data/genomeANDgff/results/prot/16cb6a/Results_Nov28"
    # prot_all_file = "/data/bioinf2023/PlantPath2023/data/genomeANDgff/results/proteinfilefinal.faa"
    # cds_all_file = "/data/bioinf2023/PlantPath2023/data/genomeANDgff/results/cdsfilefinal.fas"
    dist_matrix_kaks = geneComb(ResultsPath,arg.numberThreads,prot_all_file,cds_all_file)
    kaks_tree = dist_matrix_tree + dist_matrix_kaks

    matrix_final = append_species(kaks_tree, dict_species)
    comps = getHGT([dist_matrix_kaks, dist_matrix_tree])
    list_topology = topology(ResultsPath)
    venDiagram(comps)
    # print(comps)
    plotData(matrix_final)
    print("a")
#


