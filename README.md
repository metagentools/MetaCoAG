<p align="center">
  <img src="https://raw.githubusercontent.com/metagentools/MetaCoAG/master/MetaCoAG_Logo.png" width="500" title="MetaCoAG logo" alt="MetaCoAG logo">
</p>

# MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs

[![DOI](https://img.shields.io/badge/DOI-10.1007/978--3--031--04749--7__5-informational)](https://doi.org/10.1007/978-3-031-04749-7_5)
[![DOI](https://img.shields.io/badge/DOI-10.1089/cmb.2022.0262-green)](https://doi.org/10.1089/cmb.2022.0262)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/metacoag/badges/version.svg)](https://anaconda.org/bioconda/metacoag)
[![PyPI version](https://badge.fury.io/py/metacoag.svg)](https://badge.fury.io/py/metacoag)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/metacoag/badges/downloads.svg)](https://anaconda.org/bioconda/metacoag)

[![CI](https://github.com/metagentools/MetaCoAG/actions/workflows/testing.yml/badge.svg)](https://github.com/metagentools/MetaCoAG/actions/workflows/testing.yml)
![GitHub](https://img.shields.io/github/license/Vini2/MetaCoAG)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![CodeQL](https://github.com/metagentools/MetaCoAG/actions/workflows/codeql.yml/badge.svg)](https://github.com/metagentools/MetaCoAG/actions/workflows/codeql.yml)
[![Documentation Status](https://readthedocs.org/projects/metacoag/badge/?version=latest)](https://metacoag.readthedocs.io/en/latest/?badge=latest)

MetaCoAG is a metagenomic contig binning tool that makes use of the connectivity information found in assembly graphs, apart from the composition and coverage information. MetaCoAG makes use of single-copy marker genes along with a graph matching technique and a label propagation technique to bin contigs. MetaCoAG is tested on contigs obtained from next-generation sequencing (NGS) data. Currently, MetaCoAG supports contigs assembled using metaSPAdes and MEGAHIT, and recently we have added support for Flye assemblies (has not been tested extensively).

For detailed instructions on installation, usage and visualisation, please refer to the [documentation hosted at Read the Docs](https://metacoag.readthedocs.io/).

**NEW:** MetaCoAG is now available on bioconda at 
[https://anaconda.org/bioconda/metacoag](https://anaconda.org/bioconda/metacoag) and on PyPI at [https://pypi.org/project/metacoag/](https://pypi.org/project/metacoag/).

## Dependencies
MetaCoAG installation requires Python 3.7 or above. You will need the following python dependencies to run MetaCoAG and related support scripts. The latest tested versions of the dependencies are listed as well.
* [python](https://www.python.org/) - version 3.11.0
* [python-igraph](https://igraph.org/python/) - version 0.10.4
* [biopython](https://biopython.org/) - version 1.80
* [networkx](https://networkx.github.io/) - version 3.0
* [scipy](https://www.scipy.org/) - version 1.10.0
* [numpy](https://numpy.org/) - version 1.24.2
* [tqdm](https://github.com/tqdm/tqdm) - version 4.64.1
* [click](https://click.palletsprojects.com/) - version 8.1.3

MetaCoAG uses the following tools to scan for single-copy marker genes. These tools have been tested on the following versions.
* [FragGeneScan](https://sourceforge.net/projects/fraggenescan/) - version 1.31
* [HMMER](http://hmmer.org/) - version 3.3.2


## Installing MetaCoAG using conda

We recommend that you use [`conda`](https://docs.conda.io/en/latest/) to run MetaCoAG.

You can install MetaCoAG from [bioconda](https://anaconda.org/bioconda/metacoag).

```shell
# add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment and install metacoag
conda create -n metacoag -c bioconda metacoag

# activate metacoag environment
conda activate metacoag

# check metacoag installation
metacoag -h
```

## Setting up MetaCoAG for development

### Downloading MetaCoAG
You can clone the MetaCoAG repository to your machine.

```
git clone https://github.com/Vini2/MetaCoAG.git
```

Now go into the MetaCoAG folder using the command

```
cd MetaCoAG/
```

### Using `conda`

Once you have installed `conda`, make sure you are in the MetaCoAG folder. Now run the following commands to create a `conda` environment and activate it to run MetaCoAG.

```
conda env create -f environment.yml
conda activate metacoag
```

### Using `pip`
You can run the following command to install MetaCoAG using `pip`. Make sure you are in the MetaCoAG folder.

```
pip install .
```

**Note:** If you use pip to setup MetaCoAG for development, you will have to install [FragGeneScan](https://sourceforge.net/projects/fraggenescan/) and [HMMER](http://hmmer.org/) manually and add them to your system path.

### Test the setup

After setting up, run the following command to ensure that metacoag is working.

```
metacoag -h
```

## Example Usage

```
metacoag --assembler spades --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fasta --paths /path/to/paths_file.paths --abundance /path/to/abundance.tsv --output /path/to/output_folder
```

```
metacoag --assembler megahit --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fasta --abundance /path/to/abundance.tsv --output /path/to/output_folder
```

```
metacoag --assembler flye --graph /path/to/assembly_graph.gfa --contigs /path/to/assembly.fasta --paths /path/to/assembly_info.txt --abundance /path/to/abundance.tsv --output /path/to/output_folder
```


## Citation

MetaCoAG has been accepted at [RECOMB 2022](https://recomb2022.net/accepted-papers/) and is published in Lecture Notes in Computer Science at [https://doi.org/10.1007/978-3-031-04749-7_5](https://doi.org/10.1007/978-3-031-04749-7_5) and the journal extension is published in the Journal of Computational Biology at [https://www.liebertpub.com/doi/10.1089/cmb.2022.0262](https://www.liebertpub.com/doi/10.1089/cmb.2022.0262).

If you use MetaCoAG in your work, please cite the following publications.

```bibtex
@InProceedings{10.1007/978-3-031-04749-7_5,
  author="Mallawaarachchi, Vijini and Lin, Yu",
  editor="Pe'er, Itsik",
  title="MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs",
  booktitle="Research in Computational Molecular Biology",
  year="2022",
  publisher="Springer International Publishing",
  address="Cham",
  pages="70--85",
  abstract="Metagenomics has allowed us to obtain various genetic material from different species and gain valuable insights into microbial communities. Binning plays an important role in the early stages of metagenomic analysis pipelines. A typical pipeline in metagenomics binning is to assemble short reads into longer contigs and then bin into groups representing different species in the metagenomic sample. While existing binning tools bin metagenomic contigs, they do not make use of the assembly graphs that produce such assemblies. Here we propose MetaCoAG, a tool that utilizes assembly graphs with the composition and coverage information to bin metagenomic contigs. MetaCoAG uses single-copy marker genes to estimate the number of initial bins, assigns contigs into bins iteratively and adjusts the number of bins dynamically throughout the binning process. Experimental results on simulated and real datasets demonstrate that MetaCoAG significantly outperforms state-of-the-art binning tools, producing similar or more high-quality bins than the second-best tool. To the best of our knowledge, MetaCoAG is the first stand-alone contig-binning tool to make direct use of the assembly graph information.",
  isbn="978-3-031-04749-7"
}

@Article{doi:10.1089/cmb.2022.0262,
author = {Mallawaarachchi, Vijini and Lin, Yu},
title = {Accurate Binning of Metagenomic Contigs Using Composition, Coverage, and Assembly Graphs},
journal = {Journal of Computational Biology},
volume = {29},
number = {12},
pages = {1357--1376},
year = {2022},
doi = {10.1089/cmb.2022.0262},
note ={PMID: 36367700},
URL = {https://doi.org/10.1089/cmb.2022.0262},
eprint = {https://doi.org/10.1089/cmb.2022.0262},
abstract = { Metagenomics enables the recovery of various genetic materials from different species, thus providing valuable insights into microbial communities. Metagenomic binning group sequences belong to different organisms, which is an important step in the early stages of metagenomic analysis pipelines. The classic pipeline followed in metagenomic binning is to assemble short reads into longer contigs and then bin these resulting contigs into groups representing different taxonomic groups in the metagenomic sample. Most of the currently available binning tools are designed to bin metagenomic contigs, but they do not make use of the assembly graphs that produce such assemblies. In this study, we propose MetaCoAG, a metagenomic binning tool that uses assembly graphs with the composition and coverage information of contigs. MetaCoAG estimates the number of initial bins using single-copy marker genes, assigns contigs into bins iteratively, and adjusts the number of bins dynamically throughout the binning process. We show that MetaCoAG significantly outperforms state-of-the-art binning tools by producing similar or more high-quality bins than the second-best binning tool on both simulated and real datasets. To the best of our knowledge, MetaCoAG is the first stand-alone contig-binning tool that directly makes use of the assembly graph information along with other features of the contigs. }
}
```

## Funding

MetaCoAG is funded by an [Essential Open Source Software for Science Grant](https://chanzuckerberg.com/eoss/proposals/cogent3-python-apis-for-iq-tree-and-graphbin-via-a-plug-in-architecture/) from the Chan Zuckerberg Initiative.

<p align="left">
  <img src="https://chanzuckerberg.com/wp-content/themes/czi/img/logo.svg" width="300">
</p>
