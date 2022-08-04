<p align="center">
  <img src="MetaCoAG_Logo.png" width="500" title="MetaCoAG logo" alt="MetaCoAG logo">
</p>

# MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs

[![CI](https://github.com/metagentools/MetaCoAG/actions/workflows/testing.yml/badge.svg)](https://github.com/metagentools/MetaCoAG/actions/workflows/testing.yml)
![GitHub](https://img.shields.io/github/license/Vini2/MetaBAG)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![GitHub](https://img.shields.io/github/v/release/Vini2/MetaCoAG?include_prereleases)
[![Documentation Status](https://readthedocs.org/projects/metacoag/badge/?version=latest)](https://metacoag.readthedocs.io/en/latest/?badge=latest)

MetaCoAG is a metagenomic contig binning tool that makes use of the connectivity information found in assembly graphs, apart from the composition and coverage information. MetaCoAG makes use of single-copy marker genes along with a graph matching technique and a label propagation technique to bin contigs. MetaCoAG is tested on contigs obtained from next-generation sequencing (NGS) data. Currently, MetaCoAG supports contigs assembled using metaSPAdes and MEGAHIT, and recently we have added support for Flye assemblies (has not been tested extensively).

For detailed instructions on installation, usage and visualisation, please refer to the [documentation hosted at Read the Docs](https://metacoag.readthedocs.io/).

## Dependencies
MetaCoAG installation requires Python 3.7 or above (tested on Python 3.7.4). You will need the following python dependencies to run MetaCoAG and related support scripts. The tested versions of the dependencies are listed as well.
* [python-igraph](https://igraph.org/python/) - version 0.9.6
* [biopython](https://biopython.org/) - version 1.74
* [networkx](https://networkx.github.io/) - version 2.4
* [scipy](https://www.scipy.org/) - version 1.3.1
* [numpy](https://numpy.org/) - version 1.17.2
* [tqdm](https://github.com/tqdm/tqdm) - version 4.36.1

MetaCoAG uses the following tools to scan for single-copy marker genes. These tools are included with the following versions.
* [FragGeneScan](https://sourceforge.net/projects/fraggenescan/) - version 1.31
* [HMMER](http://hmmer.org/) - version 3.3.2


## Installing MetaCoAG

### Downloading MetaCoAG
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

### Setting up the environment
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


## Preprocessing

Firstly, you will have to assemble your set of reads into contigs. For this purpose, you can use metaSPAdes and MEGAHIT as MetaCoAG currently supports assembly graphs produced from these assemblers. Support for other assemblers will be added in future.


## Using MetaCoAG
You can see the usage options of MetaCoAG by typing `./metacoag -h` on the command line. For example,

```
usage: metacoag [-h] --assembler ASSEMBLER --graph GRAPH --contigs CONTIGS
                --abundance ABUNDANCE [--paths PATHS] --output OUTPUT
                [--hmm HMM] [--prefix PREFIX] [--min_length MIN_LENGTH]
                [--p_intra P_INTRA] [--p_inter P_INTER] [--d_limit D_LIMIT]
                [--depth DEPTH] [--mg_threshold MG_THRESHOLD]
                [--bin_mg_threshold BIN_MG_THRESHOLD]
                [--min_bin_size MIN_BIN_SIZE] [--delimiter DELIMITER]
                [--nthreads NTHREADS] [-v]

MetaCoAG is a metagenomic contig binning tool that makes use of the
connectivity information found in assembly graphs, apart from the composition
and coverage information. MetaCoAG makes use of single-copy marker genes along
with a graph matching technique and a label propagation technique to bin
contigs.

optional arguments:
  -h, --help            show this help message and exit
  --assembler ASSEMBLER
                        name of the assembler used. (Supports SPAdes, MEGAHIT
                        and Flye)
  --graph GRAPH         path to the assembly graph file
  --contigs CONTIGS     path to the contigs file
  --abundance ABUNDANCE
                        path to the abundance file
  --paths PATHS         path to the contigs.paths file
  --output OUTPUT       path to the output folder
  --hmm HMM             path to marker.hmm file. [default:
                        auxiliary/marker.hmm]
  --prefix PREFIX       prefix for the output file
  --min_length MIN_LENGTH
                        minimum length of contigs to consider for binning.
                        [default: 1000]
  --p_intra P_INTRA     minimum probability of an edge matching to assign to
                        the same bin. [default: 0.1]
  --p_inter P_INTER     maximum probability of an edge matching to create a
                        new bin. [default: 0.01]
  --d_limit D_LIMIT     distance limit for contig matching. [default: 20]
  --depth DEPTH         depth to consider for label propagation. [default: 10]
  --mg_threshold MG_THRESHOLD
                        length threshold to consider marker genes. [default:
                        0.5]
  --bin_mg_threshold BIN_MG_THRESHOLD
                        minimum fraction of marker genes that should be
                        present in a bin. [default: 0.33333]
  --min_bin_size MIN_BIN_SIZE
                        minimum size of a bin to output in base pairs.
                        [default: 200000]
  --delimiter DELIMITER
                        delimiter for output results. Supports a comma (,), a
                        semicolon (;), a tab ($'\t'), a space (" ") and a pipe
                        (|) [default: , (comma)]
  --nthreads NTHREADS   number of threads to use. [default: 8]
  -v, --version         show program's version number and exit
```


### Example Usage

```
./metacoag --assembler spades --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fasta --paths /path/to/paths_file.paths --abundance /path/to/abundance.tsv --output /path/to/output_folder
```

```
./metacoag --assembler megahit --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fasta --abundance /path/to/abundance.tsv --output /path/to/output_folder
```

```
./metacoag --assembler flye --graph /path/to/assembly_graph.gfa --contigs /path/to/assembly.fasta --paths /path/to/assembly_info.txt --abundance /path/to/abundance.tsv --output /path/to/output_folder
```


## Citation

MetaCoAG has been accepted at [RECOMB 2022](https://recomb2022.net/accepted-papers/) and is published in Lecture Notes in Computer Science at [https://doi.org/10.1007/978-3-031-04749-7_5](https://doi.org/10.1007/978-3-031-04749-7_5).

If you use MetaCoAG in your work, please cite as,

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
```
