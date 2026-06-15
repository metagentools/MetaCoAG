# Using MetaCoAG

You can see the usage options of MetaCoAG by typing `metacoag --help` on the command line. For example,

```
Usage: metacoag [OPTIONS]

  MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly
  Graphs

Options:
  --assembler [spades|megahit|megahitc|flye|custom]
                                  name of the assembler used. (Supports
                                  SPAdes, MEGAHIT and Flye)  [required]
  --graph PATH                    path to the assembly graph file  [required]
  --contigs PATH                  path to the contigs file  [required]
  --abundance PATH                path to the abundance file  [required]
  --paths PATH                    path to the contigs.paths (metaSPAdes) or
                                  assembly.info (metaFlye) file
  --output PATH                   path to the output folder  [required]
  --hmm TEXT                      path to marker.hmm file.  [default:
                                  auxiliary/marker.hmm]
  --prefix TEXT                   prefix for the output file
  --min_length INTEGER            minimum length of contigs to consider for
                                  binning.  [default: 1000]
  --p_intra FLOAT RANGE           minimum probability of an edge matching to
                                  assign to the same bin.  [default: 0.1;
                                  0<=x<=1]
  --p_inter FLOAT RANGE           maximum probability of an edge matching to
                                  create a new bin.  [default: 0.01; 0<=x<=1]
  --d_limit INTEGER               distance limit for contig matching.
                                  [default: 20]
  --depth INTEGER                 depth to consider for label propagation.
                                  [default: 10]
  --n_mg INTEGER                  total number of marker genes.  [default:
                                  108]
  --no_cut_tc                     do not use --cut_tc for hmmsearch.
  --mg_threshold FLOAT RANGE      length threshold to consider marker genes.
                                  [default: 0.5; 0<=x<=1]
  --bin_mg_threshold FLOAT RANGE  minimum fraction of marker genes that should
                                  be present in a bin.  [default: 0.33333;
                                  0<=x<=1]
  --min_bin_size INTEGER          minimum size of a bin to output in base
                                  pairs (bp).  [default: 200000]
  --delimiter [,|;|$'\t'|" "]     delimiter for output results. Supports a
                                  comma (,), a semicolon (;), a tab ($'\t'), a
                                  space (" ") and a pipe (|) .  [default: ,]
  --nthreads INTEGER              number of threads to use.  [default: 8]
  --continue                      resume from the last completed stage in the
                                  output folder.
  -v, --version                   Show the version and exit.
  --help                          Show this message and exit.
```

`min_length`, `p_intra`, `p_inter`, `d_limit`, `mg_threshold`, `bin_mg_threshold`, `min_bin_size`, `depth` and `nthreads` parameters are set by default to `1000`, `0.1`, `0.01`, `20`, `0.5`, `0.3333`, `200000`, `10` and `8` respectively. However, the user can specify them when running MetaCoAG.

If a run is interrupted, repeat the same command with `--continue`. MetaCoAG resumes after the last completed pipeline stage. The input files and result-affecting options must match the interrupted run; `--nthreads` may be changed. The checkpoint is removed automatically after a successful run.

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
metacoag --assembler spades --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fasta --paths /path/to/paths_file.paths --abundance /path/to/abundance.tsv --output /path/to/output_folder
```

```
metacoag --assembler megahit --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fasta --abundance /path/to/abundance.tsv --output /path/to/output_folder
```

```
metacoag --assembler flye --graph /path/to/assembly_graph.gfa --contigs /path/to/assembly.fasta --paths /path/to/assembly_info.txt --abundance /path/to/abundance.tsv --output /path/to/output_folder
```

# Output

The output of MetaCoAG will contain the following main files and folders.

* `contig_to_bin.tsv` containing the comma separated records of `contig id, bin number`
* `bins` containing the identified bins (FASTA file for each bin)
* `low_quality_bins` containing the identified low-quality bins, i.e., having a fraction of marker genes lower than `bin_mg_threshold` (FASTA file for each bin)
* `*.frag.faa`, `*.frag.ffn` and `*.frag.gff` files containing FragGeneScan output
* `*.hmmout` containing HMMER output
