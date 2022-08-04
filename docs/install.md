# Downloading MetaCoAG
You can download the latest release of MetaCoAG from [Releases](https://github.com/Vini2/MetaCoAG/releases) or clone the MetaCoAG repository to your machine.

```
git clone https://github.com/Vini2/MetaCoAG.git
```

If you have downloaded a release, you will have to extract the files using the following command.

```
unzip [file_name].zip
```

Now go into the MetaCoAG folder using the command

```
cd MetaCoAG/
```

# Setting up the environment
We recommend that you use [Conda](https://docs.conda.io/en/latest/) to run MetaCoAG. You can download [Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) which contains Conda.

Once you have installed Conda, make sure you are in the MetaCoAG folder. Now run the following commands to create a Conda environment and activate it to run MetaCoAG.

```
conda env create -f environment.yml
conda activate metacoag
```

Now you are ready to run MetaCoAG.

If you want to switch back to your normal environment, run the following command.

```
conda deactivate
```