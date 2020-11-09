#!/usr/bin/env python3

"""metacoag.py: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs."""

import argparse
import os
import sys
import subprocess

__author__ = "Vijini Mallawaarachchi and Yu Lin"
__copyright__ = "Copyright 2020, MetaCoAG Project"
__license__ = "GPL-3.0"
__version__ = "0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Prototype"

parser = argparse.ArgumentParser(description="""MetaCoAG is a NGS data-based metagenomic contig binning tool that makes use of the 
connectivity information found in assembly graphs, apart from the composition and coverage information. 
MetaCoAG makes use of single-copy marker genes along with a graph matching technique and a label propagation technique to bin contigs.""")

parser.add_argument("--assembler", 
                    required=True,
                    type=str,
                    help="name of the assembler used (SPAdes or MEGAHIT).")

parser.add_argument("--graph", 
                    required=True,
                    type=str,
                    help="path to the assembly graph file")

parser.add_argument("--contigs", 
                    required=True,
                    type=str,
                    help="path to the contigs file")

parser.add_argument("--abundance", 
                    required=False,
                    type=str,
                    help="path to the abundance file")

parser.add_argument("--paths", 
                    required=False,
                    type=str,
                    help="path to the contigs.paths file")

parser.add_argument("--output", 
                    required=True,
                    type=str,
                    help="path to the output folder")

parser.add_argument("--prefix", 
                    required=False,
                    type=str,
                    default='',
                    help="prefix for the output file")

parser.add_argument("--depth", 
                    required=False, 
                    type=int, 
                    default=5, 
                    help="maximum depth for the breadth-first-search. [default: 5]")

parser.add_argument("--min_length", 
                    required=False, 
                    type=int, 
                    default=1000, 
                    help="minimum length of contigs to consider for compositional probability. [default: 1000]")

parser.add_argument("--w_intra", 
                    required=False, 
                    type=int, 
                    default=2, 
                    help="maximum weight of an edge matching to assign to the same bin. [default: 2]")

parser.add_argument("--w_inter", 
                    required=False, 
                    type=int, 
                    default=80, 
                    help="minimum weight of an edge matching to create a new bin. [default: 80]")

parser.add_argument("--d_limit", 
                    required=False, 
                    type=int, 
                    default=10, 
                    help="distance limit for contig matching. [default: 10]")

parser.add_argument("--nthreads", 
                    required=False, 
                    type=int, 
                    default=8, 
                    help="number of threads to use. [default: 8]")

args = vars(parser.parse_args())


assembler = args["assembler"]
assembly_graph_file = args["graph"]
contigs = args["contigs"]
abundance = args["abundance"]
contig_paths = args["paths"]
output_path = args["output"]
prefix = args["prefix"]
depth = args["depth"]
min_length = args["min_length"]
w_intra = args["w_intra"]
w_inter = args["w_inter"]
d_limit = args["d_limit"]
nthreads = args["nthreads"]


# Validation of inputs
#---------------------------------------------------

# Check assembler name
if not (assembler.lower() == "spades"):
    print("\nPlease make sure to provide the correct assembler type.")
    print("Exiting MetaCoAG...\nBye...!\n")
    sys.exit(1)

# Check assembly graph file
if not os.path.isfile(assembly_graph_file):
    print("\nFailed to open the assembly graph file.")
    print("Exiting MetaCoAG...\nBye...!\n")
    sys.exit(1)

# Check contigs file
if not os.path.isfile(contigs):
    print("\nFailed to open the contigs file.")
    print("Exiting MetaCoAG...\nBye...!\n")
    sys.exit(1)

# Check if paths file is provided when the assembler type is SPAdes
if assembler.lower()=="spades" and contig_paths is None:
    print("\nPlease make sure to provide the path to the contigs.paths file.")
    print("Exiting MetaCoAG...\nBye...!\n")
    sys.exit(1)

# Check contigs.paths file for SPAdes
if assembler.lower()=="spades" and not os.path.isfile(contig_paths):
    print("\nFailed to open the contigs.paths file.")
    print("Exiting MetaCoAG...\nBye...!\n")
    sys.exit(1)

# Check if abundance file is provided when the assembler type is MEGAHIT
if assembler.lower() == "megahit" and abundance is None:
    print("\nPlease make sure to provide the path to the abundance file.")
    print("Exiting MetaCoAG...\nBye...!\n")
    sys.exit(1)

# Handle for missing trailing forwardslash in output folder path
if output_path[-1:] != "/":
    output_path = output_path + "/"

# Create output folder if it does not exist
if not os.path.isdir(output_path):
    subprocess.run("mkdir -p "+output_path, shell=True)

# Validate prefix
if args["prefix"] != '':
    if args["prefix"].endswith("_"):
        prefix = args["prefix"]
    else:
        prefix = args["prefix"]+"_"
else:
    prefix = ''

# Validate depth
if depth <= 0:
    print("\nPlease enter a valid number for depth")
    print("Exiting MetaCoAG...\nBye...!\n")
    sys.exit(1)

# Validate min_length
if min_length <= 0:
    print("\nPlease enter a valid number for min_length")
    print("Exiting MetaCoAG...\nBye...!\n")
    sys.exit(1)

# Validate w_intra
if w_intra <= 0:
    print("\nPlease enter a valid number for w_intra")
    print("Exiting MetaCoAG...\nBye...!\n")
    sys.exit(1)

# Validate w_inter
if w_inter <= 0:
    print("\nPlease enter a valid number for w_inter")
    print("Exiting MetaCoAG...\nBye...!\n")
    sys.exit(1)

# Validate d_limit
if d_limit <= 0:
    print("\nPlease enter a valid number for d_limit")
    print("Exiting MetaCoAG...\nBye...!\n")
    sys.exit(1)

# Validate number of threads
if nthreads <= 0:
    print("\nPlease enter a valid number for the number of threads")
    print("Exiting MetaCoAG...\nBye...!\n")
    sys.exit(1)


# Run MetaCoAG
#---------------------------------------------------

if assembler.lower()=="spades":
    cmdMetaCoAG = """python "{0}/src/metacoag_SPAdes.py" --graph "{1}" --contigs "{2}" --paths "{3}" --output "{4}" --prefix "{5}" --depth "{6}" --min_length "{7}" --w_intra "{8}" --w_inter "{9}" --d_limit "{10}" --nthreads "{11}" """.format(
        os.path.dirname(__file__), 
        assembly_graph_file,
        contigs,
        contig_paths, 
        output_path,
        prefix,
        depth,
        min_length,
        w_intra,
        w_inter,
        d_limit,
        nthreads)


os.system(cmdMetaCoAG)
