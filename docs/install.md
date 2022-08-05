# Setting up MetaCoAG

You can install MetaCoAG via [Bioconda](https://anaconda.org/bioconda/metacoag). You can download [Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) which contains Conda.

Once you have installed Conda, you can install MetaCoAG directly from the bioconda distribution using the command

```
conda install -c bioconda metacoag
```

You can also create a new conda environment and install MetaCoAG from bioconda using the following command and activate it.

```
conda create -n metacoag -c bioconda metacoag
conda activate metacoag
```

After setup, check if MetaCoAG is properly installed by typing `metacoag -h` on the command line. You should see the usage options.

If you want to switch back to your normal environment, run the following command.

```
conda deactivate
```

Now let's prepare our results to run MetaCoAG.