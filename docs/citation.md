# Citation

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