from datetime import datetime
from pathlib import Path
from itertools import combinations
from collections import defaultdict
import plotly.express as px
import argparse
from Bio import Phylo, SeqIO
import Bio
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

from Bio.SeqRecord import SeqRecord
from matplotlib_venn import venn3
import numpy as np
import shutil
import matplotlib.pyplot as plt
import ete3
from typing import Tuple, Dict

#import gffutils

# PARAAT = "ParaAT.pl  -h  %s  -a  %s  -n  %s   -p  %s  -o  %s -f axt"
# KAKS = "KaKs_Calculator  -i %s -o %s -m %s"
# ORTHOFINDER = "/data/bioinf2023/PlantPath2023/OrthoFinder/orthofinder -f %s -t %s -o %s %s"

# Docker versions
PARAAT = "ParaAT.pl -h %s -a %s -n %s -p %s -o %s -f axt -t -v"
KAKS = "KaKs -i %s -o %s -m %s"
ORTHOFINDER = "orthofinder -f %s -t %s -o %s %s"
GFFREAD = "gffread -w %s -y %s -F -S -C -g %s %s"


def arguments():
    parser = argparse.ArgumentParser(
        prog='HGT',
        description='This tool identifies potential horizontal gene transfer events between species. The input needs to be a directory containing the reference and annotation files for the species of interest. The tool will run OrthoFinder and KaKs Calculator to identify potential HGT events.', )
    parser.add_argument('-i', '--input',
                        help="In this directory should be present both reference and annotation with the same root",
                        required=True)
    parser.add_argument('-gff', '--gffread',
                        action='store_true',
                        help="Use this flag if the input files are a .fasta genome with its .gff annotation. Default is .fasta cds.")
    parser.add_argument('-OFr', '--orthofinderResults',
                        type=str,
                        help="The path to a previous OrthoFinder results folder. If not provided, the tool will run Orthofinder.")
    parser.add_argument('-nt', '--numberThreads',
                        type=int,
                        help='Number of parallel threads to use')
    parser.add_argument('-o', '--output',
                        default="output",
                        help="The name of the output file. Default is a folder named `output`.")
    parser.add_argument('-v', '--verbose',
                        action='store_true')
    parser.add_argument('-e', '--extra',
                        default='',
                        help='Eventual extra flags and attributes to pass to Orthofinder')

    args = parser.parse_args()

    return (args)

def prepare_fasta_input(arg):
    """

    :param input:
        arg: program parameters
    :return:
    """

    res_path = os.path.join(arg.input, "results")
    prot_path = os.path.join(arg.input, 'prot')
    cds_path = os.path.join(arg.input, 'cds')

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

    gene_association = {}
    irregular_proteins = []

    # Retrieve the .fasta files, containing the CDS sequences, in the input directory
    extensions = ('.fasta', '.faa', '.fa')
    fasta_files = [f for f in os.listdir(arg.input) if f.endswith(extensions)]

    assert fasta_files, 'Error: No .fasta files found in your input directory. Only .fasta file formats accepted as input.'
    assert len(fasta_files) % 2 == 0, "Error: expected 2 files per species \n\tCheck your input?"

    for filename in fasta_files:

        species_association = {}

        res_name = filename.split('.')[0]
        res_cds_file = os.path.join(cds_path, f"{res_name}_cds_mod.fas")
        res_prot_file = os.path.join(prot_path, f"{res_name}_prot_mod.faa")

        aa_seqlist = []
        cds_seqlist = []

        for cds_seq in SeqIO.parse(filename, "fasta"):

            if cds_seq.seq[0:3] != 'ATG': # Check if the sequence starts with the canonical start codon
                # if not, add it to the list of irregular proteins
                irregular_proteins.append(cds_seq)
                if arg.verbose:
                    print(f"Warning: irregular sequence, doesn't start with canonical 'ATG' start codon: {cds_seq.id}")

            # Create and assign a unique gene tag to the entry, to avoid ambiguity
            random_code = binascii.hexlify(os.urandom(10)).decode('utf8')
            random_name = f"gene_{random_code}"

            # Translate the CDS into a protein sequence
            aa_seq = cds_seq.translate(id=random_name, description='')

            # Save the gene tag and the original gene name in a dictionary
            species_association['id'] = cds_seq.id  
            species_association['species'] = res_name
            cds_seq.id = random_name
            cds_seq.description = ''

            # Append to the seqlist files
            aa_seqlist.append(cds_seq)
            cds_seqlist.append(aa_seq)

            gene_association[random_code] = species_association  # TODO check if this is correct...

        print(f'Warning: found {len(irregular_proteins)}')

        # Save the protein and CDS sequences in the respective folders
        with open(res_prot_file, "w") as Ffileaa:
            SeqIO.write(aa_seqlist, Ffileaa, "fasta")
        with open(res_cds_file, "w") as Ffilecds:
            SeqIO.write(cds_seqlist, Ffilecds, "fasta")

    gene_association_file, irregular_proteins_file = create_g_ass_and_irr_prot_files(res_path, irregular_proteins, gene_association)
    cds_all_file, prot_all_file = create_collection_file(res_path)

    return (res_path, prot_all_file, cds_all_file, gene_association_file, irregular_proteins_file)

def prepare_gff_input(arg):
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
    res_path = os.path.join(arg.input, "results")
    prot_path = os.path.join(res_path, "prot")
    cds_path = os.path.join(res_path, "cds")
    input_folder = Path(arg.input)

    # parse the gff files. this function also adds the TAG ... and modifies the gff into .gff_mod_gff

    gene_association = parse_gff(input_folder)

    # pick the modified gff files and the fasta files from the input folder.
    # (These must be the .fasta genome and the .gff annotation files)

    extensions = [".fna", ".fasta", ".gff_mod_gff", ".gff3_mod_gff", ".gtf_mod_gff"]
    filename = [file.name for file in input_folder.iterdir() if file.suffix in extensions]

    assert (len(filename) % 2) == 0, f"Error: expected 2 files per species, got {len(filename)} \n\tCheck your input?"

    try:
        if not os.path.exists(res_path):
            Path(res_path).mkdir()
        if not os.path.exists(prot_path):
            Path(prot_path).mkdir()
        if not os.path.exists(cds_path):
            Path(cds_path).mkdir()
    except FileExistsError:
        print("Warning: res folder already exist. Overwriting...")
    except Exception as e:
        print(f"Error creating directory '{res_path}': {e}")
        raise SystemExit

    # by filtering the data trough dictionary key we can run only 1 for cycle to create 
    # a nested list of file with the same name by using as a key the file name without extension

    file_dic = defaultdict(list)
    for file_name in filename:
        name = file_name.split(".")[0]
        file_dic[name].append(file_name)
    file_matched = [match for match in file_dic.values() if len(match) > 1]
    # we created a nested list with NOT sorted file that are matched

    for check in file_matched:
        assert len(check) == 2, f"Error: expected 2 files for species {check}, got {len(check)}. Check your input."

    # here the gffread run will create *.faa_mod.fasta and *.fas_mod.fasta files in the `prot` and `cds` folders
    gene_association_file, irregular_proteins_file = run_gffread(arg, file_matched)

    # Make a collection .fasta of all prot sequences and another .fasta with all cds sequences
    cds_all_file, prot_all_file = create_collection_file(res_path)

    # the output of the function is:
    #   the path to the prot folder
    #   the paths to collection files needed by KaKs (prot_all_file, cds_all_file)
    #   the gene association dictionary filepath
    #   the irregular protein filepath
    return prot_path, prot_all_file, cds_all_file, gene_association_file, irregular_proteins_file

def create_collection_file(res_path):

    prot_path = os.path.join(res_path, "prot")
    cds_path = os.path.join(res_path, "cds")

    prot_all_file = res_path + '/proteinfilefinal.faa'
    cds_all_file = res_path + '/cdsfilefinal.fas'

    # Dump all sequences into prot_all_file and
    cmd = 'cat ' + prot_path + '/*.faa_mod.fasta' + ' > ' + prot_all_file  # this .fasta file is a collection of all the aminoacid sequences
    subprocess.run(cmd, shell=True)
    cmd = 'cat ' + cds_path + '/*.fas_mod.fasta' + ' > ' + cds_all_file  # this .fasta file is a collection of all the coding DNA sequences
    subprocess.run(cmd, shell=True)

    return cds_all_file, prot_all_file

def pool_fastamod(pair: Tuple[Bio.SeqIO.SeqRecord, Bio.SeqIO.SeqRecord]) -> Tuple[
    SeqRecord, SeqRecord, Dict, bool, str]:  # parlallelized function to modify the fasta files. checks

    irregular = False
    gene_association = {}
    aa, cds = pair

    try:
        random_name = aa.description.split("HGT=")[1]
    except IndexError:
        #print(
        #    f"Error: Protein (ID = {aa.id}) has no HGT tag in the description. \nAssigning a random gene tag.")
        random_code = binascii.hexlify(os.urandom(10)).decode("utf8")
        random_name = f"gene_{random_code}"
    else:
        if ";" in random_name:
            random_name = random_name.split(";")[0]
        random_code = random_name.split("_")[-1]

    gene_association['id'] = aa.id

    aa.id = random_name
    cds.id = random_name
    aa.description = ''
    cds.description = ''

    if not aa.seq.startswith('M'):
        irregular = True  # a flag to collect all irregular proteins
    return aa, cds, gene_association, irregular, random_code

def run_gffread(arg, file_matched):
    """
    Run and parse the output of GFFREAD for each species' gff-fasta pair.

    Args:
        arg (object): The argument object containing input and other parameters.
        file_matched (list): The list of matched file names.

    Returns:
        tuple: A tuple containing the gene association dictionary and a list of irregular proteins.

    Raises:
        FileNotFoundError: If the input files are not found.

    """
    res_path = os.path.join(arg.input, "results")
    prot_path = os.path.join(res_path, "prot")
    cds_path = os.path.join(res_path, "cds")
    irregular_proteins = []
    gene_association = {}

    for x in file_matched:

        gff_f_name = str(next((f_name for f_name in x if f_name.endswith(("_mod_gff"))), None))
        genome_filename = str(next((f_name for f_name in x if f_name.endswith((".fna", ".fasta"))), None))
        res_name = str(gff_f_name.split(".")[0])

        res_cds_file = os.path.join(cds_path, f"{res_name}_cds_mod.fas")
        res_prot_file = os.path.join(prot_path, f"{res_name}_prot_mod.faa")

        fasta_genome_file = os.path.join(arg.input, genome_filename)
        gff_file = os.path.join(arg.input, gff_f_name)

        cmd = GFFREAD % (res_cds_file, res_prot_file, fasta_genome_file, gff_file)
        # GFFREAD will output the fasta with CDS and protein sequences

        result_cds = subprocess.Popen(cmd,
                                      shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE
                                      )
        out, err = result_cds.communicate()
        if arg.verbose:
            print(out)
            print(err)

        aa_seqlist = []
        cds_seqlist = []

        iterator = extract_pairs_from_files(res_cds_file,
                                            res_prot_file)  # contains all tuples (pairs of matched sequences)

        with mp.Pool(arg.numberThreads) as p:
            for x in tqdm.tqdm(p.imap_unordered(pool_fastamod, iterator), total=len(iterator),
                               desc=f"Processing '{res_prot_file}' and '{res_cds_file}' sequences..."):
                aa, cds, g_ass, irregular, random_code = x
                aa_seqlist.append(aa)
                cds_seqlist.append(cds)
                g_ass['species'] = res_name
                gene_association[random_code] = g_ass
                if irregular:
                    irregular_proteins.append(f"{gene_association[random_code]['species']}_{gene_association[random_code]['id']}")

        with open((res_prot_file + "_mod.fasta"), "w") as Ffileaa, open((res_cds_file + "_mod.fasta"), "w") as Ffilecds:
            SeqIO.write(aa_seqlist, Ffileaa, "fasta")
            SeqIO.write(cds_seqlist, Ffilecds, "fasta")

        os.remove(res_cds_file)
        os.remove(res_prot_file)

    gene_association_file, irregular_proteins_file = create_g_ass_and_irr_prot_files(res_path, irregular_proteins, gene_association)

    return gene_association_file, irregular_proteins_file

def create_g_ass_and_irr_prot_files(res_path, irregular_proteins, gene_association):

    gene_association_file = os.path.join(res_path, "gene_association.txt")

    with open(gene_association_file, "w") as GeneAssociationFile:
        for k, v in gene_association.items():
            # writing like: random_code - original id - species
            GeneAssociationFile.write(str(k) + '\t'+ str(v['id']) + '\t'+ str(v['species']) + '\n')

    irregular_proteins_file = os.path.join(res_path, "irregular_proteins.txt")

    with open(irregular_proteins_file, "w") as IrregularProteins:
        IrregularProteins.write("\n".join(irregular_proteins))

    return gene_association_file, irregular_proteins_file

def extract_pairs_from_files(res_cds_file, res_prot_file):
    iterator = []

    for aa in SeqIO.parse(res_prot_file, "fasta"):
        for cds in SeqIO.parse(res_cds_file, "fasta"):
            if aa.id == cds.id:
                pair = (aa, cds)
                iterator.append(pair)
                break
    return iterator

def parse_gff(folder):
    gene_association = {}  # gene_association is a dict used to save the mapping between the tag and the actual gene name

    gff = [".gff", ".gff3", ".gtf"]
    filename_gff = [os.path.join(folder, file.name) for file in folder.iterdir() if file.suffix in gff]

    for file in filename_gff:
        with open(file + "_mod_gff", "w") as fho:
            with open(file, "r") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    elif len(line.split("\t")) == 9 and line.split("\t")[2].startswith(("mRNA", "CDS")):
                        random_code = binascii.hexlify(os.urandom(10)).decode('utf8')
                        genes_collection = [file.rsplit(".", 1)[0].rsplit("/", 1)[1], line.replace("\t", " ").rstrip()]
                        line = line.rstrip() + ";HGT=gene_" + str(random_code)
                        fho.write(line + "\n")
                        gene_association[str(random_code)] = genes_collection
                    else:
                        fho.write(line)

    return gene_association

def run_orthofinder(prot_path, arg):
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

    # assign a path to store orthofinder results, later pass it to the orthofinder command
    orthofinder_output_folder = os.path.join(arg.input, 'results', 'OrthoFinder_results')
    # os.makedirs(orthofinder_output_folder)

    # write the linux command to run orthofinder
    runortho = ORTHOFINDER % (prot_path,
                              str(arg.numberThreads),
                              orthofinder_output_folder,
                              arg.extra)  # TODO fix eventual orthofinder extra arguments

    try:
        p1 = sp.Popen(runortho, shell=True)
        stdout, stderr = p1.communicate()
        if arg.verbose:
            print(stdout)
            print(stderr)

    except p1.returncode != 0:
        print(f"OrthoFinder failed: error code {p1.returncode}")

    # retrieve the most recent folder created in the result directory, which will contain the OrthoFinder results, and return its path
    arr = os.listdir(orthofinder_output_folder)
    orthofinder_output_folder = os.path.join(orthofinder_output_folder, arr[0])

    return orthofinder_output_folder

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
    distances = []
    genTree = Phylo.read(genetree, "newick")
    clades = list(genTree.get_terminals())
    combs = list(it.combinations(clades, 2))

    for comb in combs:

        name1 = comb[0].name
        name2 = comb[1].name

        if name1 < name2:
            name1, name2 = name2, name1

        dist = genTree.distance(name1, name2)

        type = "tree"

        name1 = name1.split("_")[-1]
        name2 = name2.split("_")[-1]

        if name1 < name2:
            name1, name2 = name2, name1

        distances.append([name1, name2, genetree.split("/")[-1].split("_")[0], float(dist), type])
        # the final output is `distances`, a list that looks like this:
        # [["gene1", "gene2", "dist", "OG", "type"]]

    return distances

def parseOrthofinder(ResultsPath: str, threads: int) -> list:
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
        for x in tqdm.tqdm(p.imap_unordered(read_tree, tree_abs_path), total=len(tree_abs_path),
                           desc="Reading GeneTrees..."):
            distances.append(x)  # appends each entry to `distances`
            pass

    distances = [item for sublist in distances for item in sublist]  # turn it into a flat list

    return distances

def kaksparallel(file: str) -> list:
    """
    Runs the KaKs calculator program in the shell.
    
    ## Args:
        file (str): The path to the .axt file to be processed.
    ## Returns:
        list_entry: A list containing the seq pair and the calculated Ks value.

    References:
    -----------
    Nei, M., & Gojobori, T. (1986). Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. Molecular biology and evolution, 3(5), 418–426. https://doi.org/10.1093/oxfordjournals.molbev.a040410
    """

    # run the 'KAKS' command in the shell
    output = file + ".kaks"
    if not os.path.exists(output) or not os.path.getsize(output) > 0:
        runkaks = KAKS % (file,  #the .axt file passed as input
                          output,  # the output file
                          "NG")  # NG is the model used for the calculation
        run = subprocess.Popen(runkaks, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = run.communicate()
        # print(out)
        print(err)

    # read the output file and return the gene pair and their distance
    list_entry = read_kaks_file(output)

    return list_entry

def read_kaks_file(kaks_filename: str) -> list:
    with open(kaks_filename, newline='') as resultskaks:
        next(resultskaks)
        for line in resultskaks:
            if line.strip():
                # take the first and fourth elements of the file
                # corresponding to the sequence name and the Ks value
                list_entry = [line.split('\t')[i] for i in [0, 3]]
    return list_entry

def parseKaKs(arg, ResultsPath, proteinfilefinal, cdsfilefinal):
    """
    A function to combine the gene pairs from the Orthofinder results and the KaKs distances.
    ## Args:
        arg: Program parameters. Duh?
        ResultsPath (str): The path to the Orthofinder results folder.
        proteinfilefinal (str): The path to the protein file.
        cdsfilefinal (str): The path to the CDS file.
    ## Returns:
        list: A list containing the gene pairs and their distances. The list will contain this data:
            `gene1` | `gene2` | `OG` | `dist` = Ks value | `type` = "kaks"
    """

    # Create the 'proc.txt' file to store the number of threads used (required by ParaAT and KaKs Calculator)
    # file_threads = os.path.join(arg.input, "proc.txt")
    file_threads = os.path.join(os.getcwd(), ResultsPath, "proc.txt")
    with open(file_threads, "w") as fh:
        fh.write(str(arg.numberThreads) + '\n')

    # Create a dictionary to match genes to orthogroups
    data = pd.read_csv(ResultsPath + "/Orthogroups/Orthogroups.tsv", sep="\t")
    data_dict = {}
    data.fillna('empty', inplace=True)

    for row in data.itertuples(index=False):
        values = []
        for i in range(1, len(data.columns)):
            key = row[0]
            values.append(row[i].split(', '))
            data_dict[key] = values

    # a dictionary that contains gene pairs matched to their orthogroup.
    # TODO questo va cambiato per evitare la run di paraAT
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

                        if gene1 < gene2:
                            gene1, gene2 = gene2, gene1

                        line = f"{gene1}\t{gene2}\n"
                        output_file.write(line)
                        # data comes in this format:
                        #   {"gene1-gene2": "OG"}
                        dict_match[gene1 + "-" + gene2] = group_name

    kaksfolder = os.path.join(os.getcwd(), 'input', 'results', "KaKs_results")

    # check if there's any .axt filepath from paraAT

    axtFiles = []  # `axtFiles will collect existing .axt alignment files
    var = []  # `var` will collect the seq pairs and their Ks value
    kaks_filepaths = []  # `kaks_filepaths` will collect existing .kaks files, from previous KaKs calculator runs

    '''if os.path.exists(kaksfolder) and os.path.getsize(kaksfolder) > 0:
        for filename in os.listdir(kaksfolder):
            if filename.endswith('.axt'):
                axtFiles.append(os.path.join(kaksfolder, filename))'''

    # else: # run ParaAT in the shell.
    # This will create .axt files from the protein and CDS files
    proteinfilefinal = os.path.join(os.getcwd(), proteinfilefinal)
    cdsfilefinal = os.path.join(os.getcwd(), cdsfilefinal)

    runparaAT = PARAAT % (file_out, proteinfilefinal, cdsfilefinal, file_threads, kaksfolder)
    # print(f"Running {runparaAT}...")
    print(file_out)
    print(proteinfilefinal)
    print(cdsfilefinal)
    print(file_threads)
    print(kaksfolder)

    run = subprocess.Popen(runparaAT, shell=True, cwd=ResultsPath, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = run.communicate()
    if arg.verbose:
        print(out)
        print(err)
    for filename in os.listdir(kaksfolder):
        if filename.endswith('.axt'):
            axtFiles.append(os.path.join(kaksfolder, filename))

    for filename in os.listdir(kaksfolder):
        if filename.endswith('.kaks'):
            kaks_filepaths.append(os.path.join(kaksfolder, filename))

    if kaks_filepaths:
        with mp.Pool(arg.numberThreads) as p:
            for x in tqdm.tqdm(p.imap_unordered(read_kaks_file, kaks_filepaths), total=len(kaks_filepaths),
                               desc="Reading .kaks files..."):
                var.append(x)  # appends each entry to `var`

    else:  # run 'kaksparallel' function on the .axt files in parallel with multiprocessing.Pool()
        with mp.Pool(arg.numberThreads) as p:
            with tqdm.tqdm(total=len(axtFiles), desc="Running KaKs Calculator...") as pbar:
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

    # create a 'distances' list with KaKs score and the gene pairs.
    #   data comes in this format:
    #   [[gene1, gene2, OG, Ks]]

    distances = []

    for entry in var:
        try:
            OG = dict_match[entry[0]]  # retrieve orthogroup by matching with the seq_pair name
        except KeyError:
            continue
        name1, name2 = entry[0].rsplit("-")  # retrieve the original gene names

        if name1 < name2:
            name1, name2 = name2, name1

        dist = entry[1]
        try:
            dist = float(dist)
        except ValueError:
            continue
        row = [name1, name2, OG, dist, "kaks"]
        distances.append(row)

    return distances

def append_species(entry_list, gene_association):
    """
    Transform an entry list (a list of rows) into a canonical pd.DataFrame.
    Add to the dataframe the information regarding the species to which the genes examined belong to
    
    Args:
    -----
    entry_list: list
        A list of lists containing the gene pairs and their distances.
    gene_association: dict
        A dictionary containing the gene names and the species to which they belong.
    
    Returns:
    --------
    pd.DataFrame
        A dataframe containing the gene pairs and their distances, with the species information added.
    """

    final = []
    for line in entry_list:

        code1 = line[0].split("_")[-1]
        code2 = line[1].split("_")[-1]

        try:
            species1 = gene_association[code1]['species']
            species2 = gene_association[code2]['species']
        except KeyError:
            print(f'Gene Code missing from association dictionary: {code1} and {code2}')
            continue
        if species1 == species2:
            continue
        if species1 < species2:
            species1, species2 = species2, species1

        row = species1 + '_vs_' + species2
        line.append(row)
        final.append(line)

    matrix = pd.DataFrame(final,
                          columns=["gene_1", "gene_2", "OG", "dist", "type", "species"])

    return matrix

def getMeanDist(comp):
    print(f"analyzing {comp.iloc[0, 4]}...")

    mean = comp['dist'].mean()
    comp["mean_dist"] = mean - comp["dist"]
    comp.sort_values(by=["mean_dist"], inplace=True, ascending=False)

    return comp[comp.mean_dist > 0]

def getHGT(matrix, gene_association):
    """
    Transform a dataframe adding the HGT bool variable; Main function to identify potential HGT events between species.
    
    Parameters:
    -----------
    `matrix` : pd.DataFrame
        A list of dataframes containing the gene pairs and their distances. The dataframes will be derived from the outputs of KaKs and OrthoFinder.

    Returns:
    --------
        `matrix` : pd.DataFrame
            A dataframe with the `HGT` column added, containing the HGT score for each gene pair.
    """

    # format the pd.DataFrame from list
    matrix2 = append_species(matrix, gene_association)

    # remove NAs
    matrix2 = matrix2[matrix2['dist'] != "NA"]  
    matrix2['dist'] = pd.to_numeric(matrix2['dist'], downcast='float')
    
    # initialize an empty column for the HGT score
    matrix2['HGT'] = None  
    sp_pairs = matrix2['species'].unique()


    # Z test implementation
    import scipy.stats as stats

    if False:
        for sp in sp_pairs:
            sp_matrix = matrix2[matrix2['species'] == sp]
            sp_matrix['dist_zscore'] = stats.zscore(matrix2['dist'])
            matrix2.loc[matrix2['species'] == sp, 'HGT'] = sp_matrix['dist_zscore'] <= -1.96

    for sp in sp_pairs: 
        # get the 5th percentile of the distance.
        threshold = matrix2[matrix2['species'] == sp]['dist'].quantile(.05)  

        # set the 'HGT' variable to True if the distance is less than the threshold
        matrix2.loc[matrix2['species'] == sp, 'HGT'] = matrix2['dist'] < threshold  

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

    dict_topology = {}  # a dictionary with OG as key and their 'HGT score' (1 or 0) as value
    list_keep = []

    for root, dir, files in os.walk(single_tree_folder):
        for file in files:
            og = file.split("_")[0]
            dict_topology[og] = 0
            file_og = os.path.join(single_tree_folder, file)
            with open(file_og, "r") as fh:
                og_single = fh.read()
                og_tree = ete3.Tree(og_single)
                for node in og_tree.traverse():
                    if node.is_leaf():
                        node.name = node.name.split("_gene")[0]
                try:
                    diff = og_tree.compare(
                        species_tree_ete)  # use the 'compare()' function from ete3 to compute the Robinson-Foulds distance
                    if diff["rf"] > 0:  # if the distance from the average species tree is greater than 0, store the 
                        list_keep.append(og)
                        dict_topology[og] = 1
                except Exception as x:
                    list_keep.append(og)
                    dict_topology[og] = 1
                    continue

    list_uniq = list(set(list_keep))  # remove duplicates from the list of HGT candidate OGs

    return list_uniq

def vennPlot(kaks_OG_list: list, tree_OG_list: list, topology_OG_list: list):
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
    fig = venn3([set(kaks_OG_list),
                 set(tree_OG_list),
                 set(topology_OG_list)],
                set_labels=('KaKs', 'Tree', 'Topology')
                )
    return fig

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
                    y="dist",
                    x="species",
                    color="type",
                    box=True,
                    hover_data=df.columns)

    return fig

def prepare_input(arg):

    # Prepare the input files, coming in a .gff format
    if arg.gffread: 
        prot_path, prot_all_file, cds_all_file, gene_association_file, irregular_proteins_file = prepare_gff_input(arg)
        gene_association = {}
        irregular_proteins = []
        with open(gene_association_file, 'r') as f:
            for i in f.readlines():
                code, name, species = i.split('\t')
                gene_association[code] = {'id': name, 'species': species}
        with open(irregular_proteins_file, 'r') as f:
            for i in f.readlines():
                irregular_proteins.append(i.strip())

    # Prepare the input files, coming in a .fasta genome format
    else: 
        prot_path, prot_all_file, cds_all_file, gene_association_file, irregular_proteins_file = prepare_fasta_input(arg)
        gene_association = {}
        irregular_proteins = []
        with open(gene_association_file, 'r') as f:
            for i in f.readlines():
                code, name, species = i.split('\t')
                gene_association[code] = {'id': name, 'species': species}
        with open(irregular_proteins_file, 'r') as f:
            for i in f.readlines():
                irregular_proteins.append(i.strip())

    return prot_all_file,cds_all_file,prot_path,gene_association,irregular_proteins

if __name__ == "__main__":

    arg = arguments()

    # Create a folder with date and time of the run to store the output
    current_date = datetime.now()
    output_folder = os.path.join(os.getcwd(), 'output', current_date.strftime("HGTResults_%d-%b-%Y_%H_%M_%S"))

    try:
        os.makedirs(output_folder)
    except FileExistsError:
        print(f"Folder '{output_folder}' already exists. Wait a second??")

    results_path = os.path.join(os.getcwd(), 'input', 'results')

    if os.path.exists(results_path) and os.listdir(results_path):

        print(f"[+] Scanning {results_path} for previous runs and results...")
        
        try:
            prot_all_file = os.path.join(results_path, "proteinfilefinal.faa")
            cds_all_file = os.path.join(results_path, "cdsfilefinal.fas")
            prot_path = os.path.join(results_path, "prot")
            cds_path = os.path.join(results_path, "cds")
            gene_association_file = os.path.join(results_path, "gene_association.txt")
            irregular_proteins_file = os.path.join(results_path, "irregular_proteins.txt")

        except FileNotFoundError:
            print('[+] Some files are missing, restarting the file preparation pipeline...')
            cmd = ['rm', '-r', results_path]
            subprocess.run(cmd, shell=True)
            prot_all_file, cds_all_file, prot_path, gene_association, irregular_proteins = prepare_input(arg)

        else:
            # create the gene_association dictionary and the irregular_proteins list from files left by the previous run
            gene_association = {}
            irregular_proteins = []
            with open(gene_association_file, 'r') as f:
                for i in f.readlines():
                    code, name, species = i.split('\t')
                    gene_association[code] = {'id': name, 'species': species}
            with open(irregular_proteins_file, 'r') as f:
                for i in f.readlines():
                    irregular_proteins.append(i.strip())

        if arg.orthofinderResults:
            orthofinder_results_path = arg.orthofinderResults

    else:

        print("[+] Reading and preparing the input files...")

        prot_all_file, cds_all_file, prot_path, gene_association, irregular_proteins = prepare_input(arg)

        print("[+] File load and preparation complete; running Orthofinder...")

    if arg.orthofinderResults:
        orthofinder_results_path = arg.orthofinderResults
    else:
        orthofinder_results_path = run_orthofinder(prot_path, arg)

    dist_matrix_tree = parseOrthofinder(orthofinder_results_path, arg.numberThreads)
    dist_matrix_tree = getHGT(dist_matrix_tree, gene_association)

    print("[+] Orthofinder scan completed; running KaKs Calculator...")

    dist_matrix_kaks = parseKaKs(arg, orthofinder_results_path, prot_all_file, cds_all_file)
    dist_matrix_kaks = getHGT(dist_matrix_kaks, gene_association)

    print("[+] KaKs Calculator run completed; checking topologies...")

    # get the list of orthogroups with significantly different topology from that of the average species tree
    list_topology = get_topology(orthofinder_results_path)

    # Create a Venn diagram of the criteria
    list_kaks = dist_matrix_kaks.loc[dist_matrix_kaks['HGT'] == True, 'OG'].to_list()
    list_tree = dist_matrix_tree.loc[dist_matrix_tree['HGT'] == True, 'OG'].to_list()
    vennpath = os.path.join(output_folder, "vennplot.png")
    fig = vennPlot(list_kaks,
                   list_tree,
                   list_topology)
    plt.savefig(vennpath)

    # Find the intersection of the three criteria
    intersection = list(set([i for i in list_kaks if i in list_tree and i in list_topology]))
    print(f'[+] Found {len(intersection)} genes that satisfy all three criteria.')
    irregular_candidates = list(set([i for i in intersection if i in irregular_proteins]))
    if irregular_candidates:
        print(f"Warning: {len(irregular_candidates)} of them are found to be irregular sequences.")


    # Go back to the original gene name from the ID introduced in the .gff file    
    match = {}
    for key, value in gene_association.items(): 
        sp_name = value['species']
        genID = value['id']
        try:
            genID = genID.split('ID=')[-1].split(";")[0]
        except:
            pass
        mvalue = '_'.join([sp_name, genID])
        match[f"gene_{key}"] = mvalue

    dist_matrix_kaks['gene_1'] = dist_matrix_kaks['gene_1'].map(match)
    dist_matrix_kaks['gene_2'] = dist_matrix_kaks['gene_2'].map(match)
    dist_matrix_tree['gene_1'] = dist_matrix_tree['gene_1'].map(match)
    dist_matrix_tree['gene_2'] = dist_matrix_tree['gene_2'].map(match)

    # Merge the two dataframes on the gene pairs
    full_matrix = pd.merge(dist_matrix_kaks,
                           dist_matrix_tree[['gene_1', 'gene_2', 'dist', 'HGT']],
                           how='inner',
                           on=['gene_1', 'gene_2'],
                           suffixes=['_kaks', '_tree'])

    # Reorder the dataframe and add the topology score
    full_matrix = full_matrix[['gene_1', 'gene_2', 'OG', 'species', 'dist_kaks', 'dist_tree', 'HGT_kaks', 'HGT_tree']]
    
    # add the topology score if the OG appears in `list_topology`
    full_matrix['HGT_topology'] = full_matrix['OG'].isin(list_topology).astype(int)
    
    # add the irregular flag if the gene is in the list of irregular proteins
    mask = (full_matrix['gene_1'].isin(irregular_proteins)) | (full_matrix['gene_2'].isin(irregular_proteins))
    full_matrix['irregular'] = mask.astype(int)

    # compute the final HGT score and sort the dataframe
    full_matrix['HGT'] = full_matrix['HGT_kaks'] + full_matrix['HGT_tree'] + full_matrix['HGT_topology']
    full_matrix.sort_values(by=['HGT'], inplace=True, ascending=False)

    print(full_matrix.head())    

    with open(os.path.join(output_folder, 'HGT_output.tsv'), 'x') as f:
        full_matrix.to_csv(f, sep='\t', index=False)
    with open(os.path.join(output_folder, 'Irregular_candidates.tsv'), 'x') as f:
        SeqIO.write(irregular_candidates, f, 'fasta')

    # Create a dataframe suitable for plotting
    import plotly
    plot_matrix = pd.concat([dist_matrix_kaks, dist_matrix_tree], axis=0)
    fig = plotData(plot_matrix)
    plotly.offline.plot(fig, filename=os.path.join(output_folder, "HGT_violin_plots.png"))