# HGT finder tool

## Workflow

The files will be processed by Orthofinder ([Emms, D.M., Kelly, S. 2019](https://doi.org/10.1186/s13059-019-1832-y)) and KaKs Calculator ([Zhang Z.](https://doi.org/10.1016/j.gpb.2021.12.002)).
Both this programs will provide values of distance between proteins of the same orthogroup (inferred by Orthofinder).
The distances that fall below the 5th quantile of the distribution will be marked for further inspection, which is made on the topology of the tree. Gene trees that diverge significantly from the inferred species tree are marked as candidates for HGT.

## Required input files:
the script needs protein sequences and their corresponding coding sequences from the investigated organism group.
The folder containing both files is passed to the program with the `-i (or --input)` argument.
It's also possible to pass a folder with previous Orthofinder results to the program to avoid multiple time-consuming runs.

Example run command:

bind-mount the directory containing your input folder to a specified input folder in the container:

`docker run -v $PWD/:/app/input cl3mente/finalhgt -i input -nt `

## Arguments

`-i or --input`: This argument is used to specify the input directory. The input directory should contain both a reference file and an annotation file with the same root.

`-OFr or --orthofinderResults`: This argument is optional and allows you to specify the OrthoFinder results file.

`-v or --verbose`: Verbose mode.

`-nt or --numberThreads`: The number of parallel threads to use.

`-o or --output`: Specifies the output directory. By default, the output will be saved in a directory named "output".
