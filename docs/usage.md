# Using MetaCoAG

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

`min_length`, `p_intra`, `p_inter`, `d_limit`, `mg_threshold`, `bin_mg_threshold`, `min_bin_size`, `depth` and `nthreads` parameters are set by default to `1000`, `0.1`, `0.01`, `20`, `0.5`, `0.3333`, `200000`, `10` and `8` respectively. However, the user can specify them when running MetaCoAG.

You can specify the delimiter for the final binning output file using the `delimiter` parameter. Enter the following values for different delimiters; 
* `,` for a comma
* `;` for a semicolon
* `$'\t'` for a tab
* `" "` for a space 
* `|` for a pipe.

# Input Format

For the metaSPAdes version, `MetaCoAG` takes in 4 files as inputs.
* Assembly graph file (in `.gfa` format)
* Contigs file (in `.fasta` format)
* Contig paths file (in `.paths` format)
* Abundance file (in `.tsv` format) with a contig in a line and its coverage in each sample separated by tabs.

For the MEGAHIT version, `MetaCoAG` takes in 3 files as inputs.
* Assembly graph file (in `.gfa` format)
* Contigs file (in `.fasta` format)
* Abundance file (in `.tsv` format) with a contig in a line and its coverage in each sample separated by tabs.

For the Flye version, `MetaCoAG` takes in 4 files as inputs.
* Assembly graph file (`assembly_graph.gfa`)
* Contigs file (`assembly.fasta`)
* Contig paths file (`assembly_info.txt`)
* Abundance file (in `.tsv` format) with a contig in a line and its coverage in each sample separated by tabs.

# Example Usage

```
./metacoag --assembler spades --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fasta --paths /path/to/paths_file.paths --abundance /path/to/abundance.tsv --output /path/to/output_folder
```

```
./metacoag --assembler megahit --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fasta --abundance /path/to/abundance.tsv --output /path/to/output_folder
```

```
./metacoag --assembler flye --graph /path/to/assembly_graph.gfa --contigs /path/to/assembly.fasta --paths /path/to/assembly_info.txt --abundance /path/to/abundance.tsv --output /path/to/output_folder
```