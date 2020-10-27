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
