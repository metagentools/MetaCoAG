# MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs

![GitHub](https://img.shields.io/github/license/Vini2/MetaBAG)
![GitHub](https://img.shields.io/github/languages/code-size/Vini2/MetaCoAG)
![GitHub](https://img.shields.io/github/v/release/Vini2/MetaCoAG?include_prereleases)

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

Firstly, you will have to assemble your set of reads into contigs. For this purpose, you can use metaSPAdes as MetaCoAG currently supports assembly graphs produced from the metaSPAdes assembler. Support for other assemblers will be added in future.

### metaSPAdes
[**SPAdes**](http://cab.spbu.ru/software/spades/) is an assembler based on the de Bruijn graph approach. [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824) is the dedicated metagenomic assembler of SPAdes. Use metaSPAdes (SPAdes in metagenomics mode) software to assemble reads into contigs. A sample command is given below.

```
spades --meta -1 Reads_1.fastq -2 Reads_2.fastq -o /path/output_folder -t 8
```

Once you have obtained the assembly output, you can run MetaCoAG.


## Using MetaCoAG
You can see the usage options of MetaCoAG by typing `./MetaCoAG -h` on the command line. For example,

```
usage: MetaCoAG [-h] --assembler ASSEMBLER --graph GRAPH --contigs CONTIGS
                [--abundance ABUNDANCE] [--paths PATHS] --output OUTPUT
                [--prefix PREFIX] [--depth DEPTH] [--min_length MIN_LENGTH]
                [--w_intra W_INTRA] [--w_inter W_INTER] [--d_limit D_LIMIT]
                [--nthreads NTHREADS]

MetaCoAG is a NGS data-based metagenomic contig binning tool that makes use of
the connectivity information found in assembly graphs, apart from the
composition and coverage information. MetaCoAG makes use of single-copy marker
genes along with a graph matching technique and a label propagation technique
to bin contigs.

optional arguments:
  -h, --help            show this help message and exit
  --assembler ASSEMBLER
                        name of the assembler used
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
  --w_intra W_INTRA     maximum weight of an edge matching to assign to the
                        same bin. [default: 2]
  --w_inter W_INTER     minimum weight of an edge matching to create a new
                        bin. [default: 80]
  --d_limit D_LIMIT     distance limit for contig matching. [default: 10]
  --nthreads NTHREADS   number of threads to use. [default: 8]
```

`min_length`, `w_intra`, `w_inter`, `d_limit`, and `nthreads` parameters are set by default to `1000`, `2`, `80`, `10`, and `8` respectively. However, the user can specify them when running MetaCoAG.

## Input Format

For the metaSPAdes version, `metacoag` takes in 3 files as inputs (required).
* Assembly graph file (in `.gfa` format)
* Contigs file (in `.fasta` format)
* Paths of contigs (in `.paths` format)

## Example Usage

```
./MetaCoAG --assembler spades --graph /path/to/graph_file.gfa --paths /path/to/paths_file.paths --output /path/to/output_folder
```

## Test Data

Three metaSPAdes datasets used to test MetaCoAG can be found in the `test_data` folder. The test data for each of the datasets include the following files.
* Contigs file
* Assembly graph file
* Contigs file
* Paths file for the assembly graph
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
python evaluate.py --binned /path/to/binning_result.csv --groundtruth /path/to/ground_truth.csv
```

## References

[1] Albertsen,  M.,  Hugenholtz,  P.,  Skarshewski,  A.,  Nielsen,  K.L.,  Tyson,  G.W.,Nielsen,  P.H.:  Genome  sequences  of  rare,  uncultured  bacteria  obtained  by  dif-ferential coverage binning of multiple metagenomes. Nature Biotechnology 31(6),533–538 (Jun 2013)

[2] Alneberg,  J.,  Bjarnason,  B.S.,  de  Bruijn,  I.,  Schirmer,  M.,  Quick,  J.,  Ijaz,  U.Z.,Lahti, L., Loman, N.J., Andersson, A.F., Quince, C.: Binning metagenomic contigsby coverage and composition. Nature Methods 11, 1144–1146 (Sep 2014)

[3] Bankevich, A., Nurk, S., Antipov, D., Gurevich, A.A., Dvorkin, M., Kulikov, A.S.,Lesin, V.M., Nikolenko, S.I., Pham, S., Prjibelski, A.D., Pyshkin, A.V., Sirotkin,A.V.,  Vyahhi,  N.,  Tesler,  G.,  Alekseyev,  M.A.,  Pevzner,  P.A.:  SPAdes:  A  NewGenome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Jour-nal of Computational Biology 19(5), 455–477 (2012), pMID: 22506599

[4] Barnum, T.P., Figueroa, I.A., Carlstr ̈om, C.I., Lucas, L.N., Engelbrektson, A.L.,Coates, J.D.: Genome-resolved metagenomics identifies genetic mobility, metabolicinteractions, and unexpected diversity in perchlorate-reducing communities. The ISME Journal 12(6), 1568–1581 (2018)

[5] Kang, D., Li, F., Kirton, E.S., Thomas, A., Egan, R.S., An, H., Wang, Z.: MetaBAT2:  an  adaptive  binning  algorithm  for  robust  and  efficient  genome  reconstructionfrom metagenome assemblies. PeerJ7, e27522v1 (Feb 2019)

[6] Karp, R.M.: An algorithm to solve the m×n assignment problem in expectedtime O(mn log n). Networks10(2), 143–152 (1980)

[7] Kolmogorov, M., Yuan, J., Lin, Y., Pevzner, P.A.: Assembly of long, error-pronereads using repeat graphs. Nature biotechnology 37(5), 540–546 (2019)

[8] Laczny, C.C., Kiefer, C., Galata, V., Fehlmann, T., Backes, C., Keller, A.: BusyBeeWeb: metagenomic data analysis by bootstrapped supervised binning and annota-tion. Nucleic Acids Research 45(W1), W171–W179 (05 2017)

[9] Lin, H.H., Liao, Y.C.: Accurate binning of metagenomic contigs via automated clustering sequences using information of genomic signatures and marker genes. Scientific Reports 6(1), 24175 (Apr 2016)

[10] Mallawaarachchi, V., Wickramarachchi, A., Lin, Y.: GraphBin: Refined binning ofmetagenomic contigs using assembly graphs. Bioinformatics 36(11), 3307-3313, (03 2020)

[11] Mallawaarachchi, V.G., Wickramarachchi, A.S., Lin, Y.: GraphBin2: Refined andOverlapped Binning of Metagenomic Contigs Using Assembly Graphs. In: Kings-ford, C., Pisanti, N. (eds.) 20th International Workshop on Algorithms in Bioinfor-matics (WABI 2020). Leibniz International Proceedings in Informatics (LIPIcs),vol. 172, pp. 8:1–8:21. Schloss Dagstuhl–Leibniz-Zentrum f ̈ur Informatik, Dagstuhl,Germany (2020)

[12] Myers,  E.W.:  The  fragment  assembly  string  graph.  Bioinformatics 21(suppl2), ii79–ii85 (09 2005)

[13] [NetworkX: networkx.algorithms.bipartite.matching.minimumweightfullmatching - NetworkX 2.5 documentation](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.bipartite.matching.minimum_weight_full_matching.html) (2020)

[14] Nurk, S., Meleshko, D., Korobeynikov, A., Pevzner, P.A.: metaSPAdes: a newversatile metagenomic assembler. Genome Research 27(5), 824–834 (2017)

[15] Strous, M., Kraft, B., Bisdorf, R., Tegetmeyer, H.: The Binning of MetagenomicContigs for Microbial Physiology of Mixed Cultures. Frontiers in Microbiology 3,410 (2012)

[16] Wu, Y.W., Simmons, B.A., Singer, S.W.: MaxBin 2.0: an automated binning al-gorithm to recover genomes from multiple metagenomic datasets. Bioinformatics 32(4), 605–607 (Oct 2015)

[17] Yu, G., Jiang, Y., Wang, J., Zhang, H., Luo, H.: BMC3C: binning metagenomiccontigs using codon usage, sequence composition and read coverage. Bioinformatics 34(24), 4172–4179 (06 2018)
