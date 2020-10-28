# MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs

![GitHub](https://img.shields.io/github/license/Vini2/MetaBAG) 
![GitHub top language](https://img.shields.io/github/languages/top/Vini2/MetaBAG)

MetaCoAG is a NGS data-based metagenomic contig binning tool that makes use of theconnectivity information found in assembly graphs, apart from the composition and coverage information. MetaCoAG makes use of single-copy marker genes along with a graph matching technique and a label propagation technique to bin contigs.

## Dependencies
MetaCoAG installation requires Python 3 (tested on Python 3.7.4). You will need the following python dependencies to run MetaCoAG and related support scripts. The tested versions of the dependencies are listed as well.
* [python-igraph](https://igraph.org/python/) - version 0.7.1
* [biopython](https://biopython.org/) - version 1.74
* [networkx](https://networkx.github.io/) - version 2.4
* [scipy](https://www.scipy.org/) - version 1.3.1
* [numpy](https://numpy.org/) - version 1.17.2
* [tqdm](https://github.com/tqdm/tqdm) - version 4.36.1
* [tabulate](https://pypi.org/project/tabulate/) - version 0.8.7 (for [`evaluate.py`](https://github.com/Vini2/MetaCoAG/blob/master/evaluation_scripts/evaluate.py))

MetaCoAG uses the following tools to scan for single-copy marker genes. These tools are included with the following versions.
* [FragGeneScan](https://sourceforge.net/projects/fraggenescan/) - version 1.31
* [HMMER](http://hmmer.org/) - version 3.3


### Downloading MetaCoAG
You can download the latest release of MetaCoAG from [Releases](https://github.com/Vini2/MetaCoAG/releases) or clone the MetaCoAG repository to your machine.

```
git clone https://github.com/Vini2/MetaCoAG.git
```

If you have downloaded a release, you will have to extract the files using the following command.

```
unzip [file_name].zip
```

Now go in to the MetaCoAG folder using the command

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

Firstly, you will have to assemble your set of reads into contigs. For this purpose, you can use metaSPAdes or MEGAHIT as MetaCoAG currently supports assembly graphs produced from these two assemblers.

### metaSPAdes
[**SPAdes**](http://cab.spbu.ru/software/spades/) is an assembler based on the de Bruijn graph approach. [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824) is the dedicated metagenomic assembler of SPAdes. Use metaSPAdes (SPAdes in metagenomics mode) software to assemble reads into contigs. A sample command is given below.

```
spades --meta -1 Reads_1.fastq -2 Reads_2.fastq -o /path/output_folder -t 16
```

### MEGAHIT
[MEGAHIT](https://github.com/voutcn/megahit) is an assembler based on the de Bruijn graph approach. Use MEGAHIT software to assemble reads into contigs. A sample command is given below.

Once you have obtained the assembly output, you can run MetaCoAG.


## Using MetaCoAG
You can see the usage options of MetaCoAG by typing `./metacoag.py -h` on the command line. For example,

```
usage: metacoag.py [-h] --assembler ASSEMBLER --graph GRAPH --contigs CONTIGS
                   [--abundance ABUNDANCE] [--paths PATHS] --output OUTPUT
                   [--prefix PREFIX] [--depth DEPTH] [--min_length MIN_LENGTH]
                   [--alpha_intra ALPHA_INTRA] [--alpha_inter ALPHA_INTER]
                   [--dist_intra DIST_INTRA] [--dist_inter DIST_INTER]
                   [--nthreads NTHREADS]

MetaCoAG is a NGS data-based metagenomic contig binning tool that makes use of
the connectivity information found in assembly graphs, apart from the
composition and coverage information. MetaCoAG makes use of single-copy marker
genes along with a graph matching technique and a label propagation technique
to bin contigs.

optional arguments:
  -h, --help            show this help message and exit
  --assembler ASSEMBLER
                        name of the assembler used (SPAdes or MEGAHIT).
  --graph GRAPH         path to the assembly graph file
  --contigs CONTIGS     path to the contigs file
  --abundance ABUNDANCE
                        path to the abundance file
  --paths PATHS         path to the contigs.paths file
  --output OUTPUT       path to the output folder
  --prefix PREFIX       prefix for the output file
  --depth DEPTH         maximum depth for the breadth-first-search. [default:
                        5]
  --min_length MIN_LENGTH
                        minimum length of contigs to consider for
                        compositional probability. [default: 1000]
  --alpha_intra ALPHA_INTRA
                        maximum weight of an edge matching to assign to the
                        same bin. [default: 2]
  --alpha_inter ALPHA_INTER
                        minimum weight of an edge matching to create a new
                        bin. [default: 80]
  --dist_intra DIST_INTRA
                        maximum distance of a contig matched to assign to the
                        same bin. [default: 10]
  --dist_inter DIST_INTER
                        minimum distance of a contig matched to create a new
                        bin. [default: 10]
  --nthreads NTHREADS   number of threads to use. [default: 8]
```

`min_length`, `alpha_intra`, `alpha_inter`, `dist_intra`, `dist_inter` and `nthreads` parameters are set by default to `1000`, `2`, `80`, `10`, `10` and `8` respectively. However, the user can specify them when running MetaCoAG.

## Input Format

For the metaSPAdes version, `metacoag` takes in 3 files as inputs (required).
* Assembly graph file (in `.gfa` format)
* Contigs file (in `.fasta` format)
* Paths of contigs (in `.paths` format)

For the MEGAHIT version, `metacoag` takes in 3 files as inputs (required).
* Assembly graph file (in `.gfa` format. To convert fastg to gfa refer [here](https://github.com/Vini2/GraphBin/blob/master/support/README.md#fastg2gfa))
* Contigs file (in `.fa` format)
* Abundance file (tab separated file with contig ID and coverage in each line)


## Example Usage

```
./metacoag.py --assembler spades --graph /path/to/graph_file.gfa --paths /path/to/paths_file.paths --output /path/to/output_folder
```
```
./metacoag.py --assembler megahit --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fa --abundance /path/to/abundance.tsv --output /path/to/output_folder
```

## Test Data

The data used to test MetaCoAG can be found in the `test_data` folder. The test data for each of the datasets include the following files.
* Contigs file
* Assembly graph file
* Contigs file
* Paths file for the assembly graph (for the datasets assembled using metaSPAdes)
* Abundance file (for the datasets assembled using MEGAHIT)
* Initial binning result from [MaxBin 2.0](https://sourceforge.net/projects/maxbin2/)
* Initial binning result from [MetaWatt](https://sourceforge.net/p/metawatt/wiki/Home/)
* Initial binning result from [CONCOCT](https://concoct.readthedocs.io/en/latest/)
* Initial binning result from [MetaBAT 2](https://bitbucket.org/berkeleylab/metabat/src/master/)
* Initial binning result from [SolidBin](https://github.com/sufforest/SolidBin)
* Initial binning result from [BusyBee Web](https://ccb-microbe.cs.uni-saarland.de/busybee/) (Not available for metaSPAdes assemblies)
* Ground truth labelling of contigs

You can try running MetaCoAG using these test data files.

## Support Scripts

MetaCoAG provides the evaluation script [`evaluate.py`](https://github.com/Vini2/MetaCoAG/blob/master/evaluation_scripts/evaluate.py) to evaluate binning results against a known ground truth.

```
python evaluate.py --binned /path/to/binning_result.csv --groundtruth /path/to/graound_truth.csv
```
