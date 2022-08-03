# Introduction

## MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs

MetaCoAG is a metagenomic contig binning tool that makes use of the connectivity information found in assembly graphs, apart from the composition and coverage information. MetaCoAG makes use of single-copy marker genes along with a graph matching technique and a label propagation technique to bin contigs. MetaCoAG is tested on contigs obtained from next-generation sequencing (NGS) data. Currently, MetaCoAG supports contigs assembled using metaSPAdes and MEGAHIT, and recently we have added support for Flye assemblies (has not been tested extensively).
