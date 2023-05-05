# Preprocessing

## Assembly

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
**Note:** Currently, MetaCoAG supports GFA file format for the assembly graph file. The MEGAHIT toolkit will produce a FASTG file which you can convert to GFA format using [fastg2gfa](https://github.com/lh3/gfa1/blob/master/misc/fastg2gfa.c).

```
fastg2gfa final.fastg > final.gfa
```
Support for FASTG files will be added in the near future.

### Flye
[**Flye**](https://github.com/fenderglass/Flye) is a long-read assembler based on the de Bruijn graph approach. **metaFlye** is the metagenomic version of Flye. Use metaFlye to assemble reads into contigs. A sample command is given below.

```
flye --meta --pacbio-raw Reads.fastq --out-dir /path/output_folder --threads 8
```

## How to get the abundance.tsv file

You can use [CoverM](https://github.com/wwood/CoverM) to get the coverage of contigs. You can run the following commands to get the `abundance.tsv` file.

```
coverm contig -1 reads_1.fastq -2 reads_2.fastq -r contigs.fasta -o abundance.tsv -t 8
sed -i '1d' abundance.tsv	# remove the header of the file
```

You can use the -c (or --coupled) option of CoverM if you have multiple samples. Please refer the [CoverM contig documentation](https://wwood.github.io/CoverM/coverm-contig.html) for further details.

The resulting `abundance.tsv` file can be directly used in MetaCoAG.

Once you have obtained the assembly output and the `abundance.tsv` file, you can run MetaCoAG.