Example run command:

bind-mount the directory containing your dataset to a specified input folder in the container:

`docker run --rm -it --user $(id -u):$(id -g) -v $PWD/:/app/input prova1 bash`