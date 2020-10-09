# MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs

![GitHub](https://img.shields.io/github/license/Vini2/MetaBAG) 
![GitHub top language](https://img.shields.io/github/languages/top/Vini2/MetaBAG)

MetaCoAG is a NGS data-based metagenomic contig binning tool that makes use of theconnectivity information found in assembly graphs, apart from the composition and coverage information. MetaCoAG makes use of single-copy marker genes along with a graph matching technique and a label propagation technique to bin contigs.

## Dependencies
MetaCoAG installation requires python 3.6 or above. You will need the following dependencies to run MetaCoAG and related support scripts.
* [python-igraph](https://igraph.org/python/)
* [biopython](https://biopython.org/)
* [networkx](https://networkx.github.io/)
* [scipy](https://www.scipy.org/)
* [numpy](https://numpy.org/)
* [tqdm](https://github.com/tqdm/tqdm)
