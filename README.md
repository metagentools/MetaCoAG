<p align="center">
  <img src="MetaCoAG_Logo.png" width="500" title="MetaCoAG logo" alt="MetaCoAG logo">
</p>

# MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs

![GitHub](https://img.shields.io/github/license/Vini2/MetaBAG)
![GitHub](https://img.shields.io/github/languages/code-size/Vini2/MetaCoAG)
![GitHub](https://img.shields.io/github/v/release/Vini2/MetaCoAG?include_prereleases)

MetaCoAG is a metagenomic contig binning tool that makes use of the connectivity information found in assembly graphs, apart from the composition and coverage information. MetaCoAG makes use of single-copy marker genes along with a graph matching technique and a label propagation technique to bin contigs. MetaCoAG is tested on contigs obtained from next-generation sequencing (NGS) data. Currently MetaCoAG supports contigs assembled using metaSPAdes and MEGAHIT.

## Dependencies
MetaCoAG installation requires Python 3.7 (tested on Python 3.7.4). You will need the following python dependencies to run MetaCoAG and related support scripts. The tested versions of the dependencies are listed as well.
* [python-igraph](https://igraph.org/python/) - version 0.9.6
* [biopython](https://biopython.org/) - version 1.74
* [networkx](https://networkx.github.io/) - version 2.4
* [scipy](https://www.scipy.org/) - version 1.3.1
* [numpy](https://numpy.org/) - version 1.17.2
* [tqdm](https://github.com/tqdm/tqdm) - version 4.36.1

MetaCoAG uses the following tools to scan for single-copy marker genes. These tools are included with the following versions.
* [FragGeneScan](https://sourceforge.net/projects/fraggenescan/) - version 1.31
* [HMMER](http://hmmer.org/) - version 3.3.2


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

Firstly, you will have to assemble your set of reads into contigs. For this purpose, you can use metaSPAdes and MEGAHIT as MetaCoAG currently supports assembly graphs produced from these assemblers. Support for other assemblers will be added in future.

### metaSPAdes
[**SPAdes**](http://cab.spbu.ru/software/spades/) is an assembler based on the de Bruijn graph approach. [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824) is the dedicated metagenomic assembler of SPAdes. Use metaSPAdes (SPAdes in metagenomics mode) software to assemble reads into contigs. A sample command is given below.

```
spades --meta -1 Reads_1.fastq -2 Reads_2.fastq -o /path/output_folder -t 8
```

### MEGAHIT
[**MEGAHIT**](https://github.com/voutcn/megahit) is an assembler based on the de Bruijn graph approach. Use MEGAHIT software to assemble reads into contigs. A sample command is given below.

```
megahit -1 Reads_1.fastq -2 Reads_2.fastq --k-min 21 --k-max 77 -o /path/output_folder -t 8
```
**Note:** Currently, MetaCoAG supports GFA file format for the assembly graph file. The MEGAHIT toolkit will result in a FASTG file which you can convert to GFA format using [fastg2gfa](https://github.com/lh3/gfa1/blob/master/misc/fastg2gfa.c).

```
fastg2gfa final.fastg > final.gfa
```
Support for FASTG files will be added in the near future.

Once you have obtained the assembly output, you can run MetaCoAG.


## Using MetaCoAG
You can see the usage options of MetaCoAG by typing `./MetaCoAG -h` on the command line. For example,

```
usage: MetaCoAG [-h] --assembler ASSEMBLER --graph GRAPH --contigs CONTIGS
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

`min_length`, `p_intra`, `p_inter`, `d_limit`, `mg_threshold`, `bin_mg_threshold`, `min_bin_size`, `depth` and `nthreads` parameters are set by default to `1000`, `0.1`, `0.01`, `20`, `0.5`, `0.3333`, `200000`, `10` and `8` respectively. However, the user can specify them when running MetaCoAG.

You can specify the delimiter for the final binning output file using the `delimiter` paramter. Enter the following values for different delimiters; 
* `,` for a comma
* `;` for a semicolon
* `$'\t'` for a tab
* `" "` for a space 
* `|` for a pipe.

## Input Format

For the metaSPAdes version, `MetaCoAG` takes in 4 files as inputs.
* Assembly graph file (in `.gfa` format)
* Contigs file (in `.fasta` format)
* Contig paths file (in `.paths` format)
* Abundance file (in `.tsv` format) with a contig in a line and its coverage in each sample.

For the MEGAHIT version, `MetaCoAG` takes in 3 files as inputs.
* Assembly graph file (in `.gfa` format)
* Contigs file (in `.fasta` format)
* Abundance file (in `.tsv` format) with a contig in a line and its coverage in each sample.

## Example Usage

```
./MetaCoAG --assembler spades --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fasta --paths /path/to/paths_file.paths --abundance /path/to/abundance.tsv --output /path/to/output_folder
```

## References

[1] Albertsen,  M.,  Hugenholtz,  P.,  Skarshewski,  A.,  Nielsen,  K.L.,  Tyson,  G.W.,Nielsen,  P.H.:  Genome  sequences  of  rare,  uncultured  bacteria  obtained  by  dif-ferential coverage binning of multiple metagenomes. Nature Biotechnology 31(6),533–538 (Jun 2013)

[2] Alneberg,  J.,  Bjarnason,  B.S.,  de  Bruijn,  I.,  Schirmer,  M.,  Quick,  J.,  Ijaz,  U.Z.,Lahti, L., Loman, N.J., Andersson, A.F., Quince, C.: Binning metagenomic contigsby coverage and composition. Nature Methods 11, 1144–1146 (Sep 2014)

[3] Bankevich, A., Nurk, S., Antipov, D., Gurevich, A.A., Dvorkin, M., Kulikov, A.S.,Lesin, V.M., Nikolenko, S.I., Pham, S., Prjibelski, A.D., Pyshkin, A.V., Sirotkin,A.V.,  Vyahhi,  N.,  Tesler,  G.,  Alekseyev,  M.A.,  Pevzner,  P.A.:  SPAdes:  A  NewGenome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Jour-nal of Computational Biology 19(5), 455–477 (2012), pMID: 22506599

[4] Barnum, T.P., Figueroa, I.A., Carlstr ̈om, C.I., Lucas, L.N., Engelbrektson, A.L.,Coates, J.D.: Genome-resolved metagenomics identifies genetic mobility, metabolicinteractions, and unexpected diversity in perchlorate-reducing communities. The ISME Journal 12(6), 1568–1581 (2018)

[5] Kang, D., Li, F., Kirton, E.S., Thomas, A., Egan, R.S., An, H., Wang, Z.: MetaBAT2:  an  adaptive  binning  algorithm  for  robust  and  efficient  genome  reconstructionfrom metagenome assemblies. PeerJ7, e27522v1 (Feb 2019)

[6] Karp, R.M.: An algorithm to solve the m×n assignment problem in expectedtime O(mn log n). Networks10(2), 143–152 (1980)

[7] Kolmogorov, M., Yuan, J., Lin, Y., Pevzner, P.A.: Assembly of long, error-pronereads using repeat graphs. Nature biotechnology 37(5), 540–546 (2019)

[9] Lin, H.H., Liao, Y.C.: Accurate binning of metagenomic contigs via automated clustering sequences using information of genomic signatures and marker genes. Scientific Reports 6(1), 24175 (Apr 2016)

[10] Mallawaarachchi, V., Wickramarachchi, A., Lin, Y.: GraphBin: Refined binning ofmetagenomic contigs using assembly graphs. Bioinformatics 36(11), 3307-3313, (03 2020)

[11] Mallawaarachchi, V.G., Wickramarachchi, A.S., Lin, Y.: GraphBin2: Refined andOverlapped Binning of Metagenomic Contigs Using Assembly Graphs. In: Kings-ford, C., Pisanti, N. (eds.) 20th International Workshop on Algorithms in Bioinfor-matics (WABI 2020). Leibniz International Proceedings in Informatics (LIPIcs),vol. 172, pp. 8:1–8:21. Schloss Dagstuhl–Leibniz-Zentrum f ̈ur Informatik, Dagstuhl,Germany (2020)

[12] Myers,  E.W.:  The  fragment  assembly  string  graph.  Bioinformatics 21(suppl2), ii79–ii85 (09 2005)

[13] [NetworkX: networkx.algorithms.bipartite.matching.minimumweightfullmatching - NetworkX 2.5 documentation](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.bipartite.matching.minimum_weight_full_matching.html) (2020)

[14] Nissen, J. N., et al. Improved metagenome binning and assembly using deep variational autoencoders. Nature Biotechnology 39, 555–560 (2021)

[15] Nurk, S., Meleshko, D., Korobeynikov, A., Pevzner, P.A.: metaSPAdes: a newversatile metagenomic assembler. Genome Research 27(5), 824–834 (2017)

[16] Wu, Y.W., Simmons, B.A., Singer, S.W.: MaxBin 2.0: an automated binning al-gorithm to recover genomes from multiple metagenomic datasets. Bioinformatics 32(4), 605–607 (Oct 2015)

## Citation

MetaCoAG has been accepted at [RECOMB 2022](https://recomb2022.net/accepted-papers/) and the preprint is available at bioRxiv ([https://doi.org/10.1101/2021.09.10.459728](https://doi.org/10.1101/2021.09.10.459728)).

If you use MetaCoAG in your work, please cite as,

```bibtex
@article {Mallawaarachchi2021.09.10.459728,
	author = {Mallawaarachchi, Vijini and Lin, Yu},
	title = {MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs},
	elocation-id = {2021.09.10.459728},
	year = {2021},
	doi = {10.1101/2021.09.10.459728},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Metagenomics binning has allowed us to study and characterize various genetic material of different species and gain insights into microbial communities. While existing binning tools bin metagenomics de novo assemblies, they do not make use of the assembly graphs that produce such assemblies. Here we propose MetaCoAG, a tool that utilizes assembly graphs with the composition and coverage information to bin metagenomic contigs. MetaCoAG uses single-copy marker genes to estimate the number of initial bins, assigns contigs into bins iteratively and adjusts the number of bins dynamically throughout the binning process. Experimental results on simulated and real datasets demonstrate that MetaCoAG significantly outperforms state-of-the-art binning tools, producing more high-quality bins than the second-best tool, with an average median F1-score of 88.40\%. To the best of our knowledge, MetaCoAG is the first stand-alone binning tool to make direct use of the assembly graph information. MetaCoAG is available at https://github.com/Vini2/MetaCoAG.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2021/09/11/2021.09.10.459728},
	eprint = {https://www.biorxiv.org/content/early/2021/09/11/2021.09.10.459728.full.pdf},
	journal = {bioRxiv}
}
```
