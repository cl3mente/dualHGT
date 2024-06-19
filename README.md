# HGT finder tool
[![Static Badge](https://img.shields.io/badge/DockerHub-dualHGT-gray?style=flat&logo=docker&logoColor=white&labelColor=blue&link=https%3A%2F%2Fhub.docker.com%2Frepository%2Fdocker%2Fcl3mente%2Fdualhgt%2Fgeneral)](https://hub.docker.com/repository/docker/cl3mente/dualhgt/general)
[![GitHub commit activity](https://img.shields.io/github/commit-activity/t/cl3mente/dualHGT?logo=github&logoColor=black&labelColor=white&color=gray)](https://github.com/cl3mente/dualHGT/)

## Workflow

The files will be processed by Orthofinder ([Emms, D.M., Kelly, S. 2019](https://doi.org/10.1186/s13059-019-1832-y)) and KaKs Calculator ([Zhang Z.](https://doi.org/10.1016/j.gpb.2021.12.002)).
Both this programs will provide values of distance between proteins of the same orthogroup (inferred by Orthofinder).
The distances that fall below the 5th quantile of the distribution will be marked for further inspection, which is made on the topology of the tree. Gene trees that diverge significantly from the inferred species tree are marked as candidates for HGT.

## Usage:
### DockerHub pull:
Run this in your terminal to pull a pre-made image from DockerHub:

`docker pull cl3mente/dualHGT:latest`

After the pull is complete, you will have a working image that you can use to run dualHGT on a container - specify your local input folder binding it to an _'input'_ volume in the container:

`docker run -v your-input-folder/:/app/input cl3mente/dualhgt:latest`

Alternatively, run this command from the folder itself:

`docker run -v $PWD/:/app/input cl3mente/dualhgt:latest`

This folder will be the channel for Docker to communicate with your local machine. This is where the results will be written, in the _'output'_ subfolder.
Once you run the container, run the dualHGT.py script with this command and the other options that you might want to customize:

`python dualHGT.py -i input/ [...]`

### Additional arguments

`-i or --input`: This argument is used to specify the input directory. 

`-gff`: A flag specifying whether the input directory contains:
 - FLAG ON: a reference genome .fasta and an annotation .gff with the same root for each species selected (for this)
 - FLAG OFF (default): multifasta files with protein sequences of the species selected

`-OFr or --orthofinderResults`: This argument is optional and allows you to specify the OrthoFinder results file.

`-v or --verbose`: Verbose mode.

`-nt or --numberThreads`: The number of threads to use for the analysis.

the script needs protein sequences and their corresponding coding sequences from the investigated organism group.
The folder containing both files is passed to the program with the `-i (or --input)` argument.
It's also possible to pass a folder with previous Orthofinder results to the program to avoid multiple time-consuming runs.
