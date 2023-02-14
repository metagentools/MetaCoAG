#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys


def get_args(version):
    parser = argparse.ArgumentParser(
        description="""MetaCoAG is a metagenomic contig binning tool that makes use of the 
    connectivity information found in assembly graphs, apart from the composition and coverage information. 
    MetaCoAG makes use of single-copy marker genes along with a graph matching technique and a label propagation technique to bin contigs."""
    )

    parser.add_argument(
        "--assembler",
        required=True,
        type=str,
        help="name of the assembler used. (Supports SPAdes, MEGAHIT and Flye)",
    )

    parser.add_argument(
        "--graph", required=True, type=str, help="path to the assembly graph file"
    )

    parser.add_argument(
        "--contigs", required=True, type=str, help="path to the contigs file"
    )

    parser.add_argument(
        "--abundance", required=True, type=str, help="path to the abundance file"
    )

    parser.add_argument(
        "--paths",
        required=False,
        type=str,
        help="path to the contigs.paths (metaSPAdes) or assembly.info (metaFlye) file",
    )

    parser.add_argument(
        "--output", required=True, type=str, help="path to the output folder"
    )

    parser.add_argument(
        "--hmm",
        required=False,
        type=str,
        default="",
        help="path to marker.hmm file. [default: auxiliary/marker.hmm]",
    )

    parser.add_argument(
        "--prefix",
        required=False,
        type=str,
        default="",
        help="prefix for the output file",
    )

    parser.add_argument(
        "--min_length",
        required=False,
        type=int,
        default=1000,
        help="minimum length of contigs to consider for binning. [default: 1000]",
    )

    parser.add_argument(
        "--p_intra",
        required=False,
        type=float,
        default=0.1,
        help="minimum probability of an edge matching to assign to the same bin. [default: 0.1]",
    )

    parser.add_argument(
        "--p_inter",
        required=False,
        type=float,
        default=0.01,
        help="maximum probability of an edge matching to create a new bin. [default: 0.01]",
    )

    parser.add_argument(
        "--d_limit",
        required=False,
        type=int,
        default=20,
        help="distance limit for contig matching. [default: 20]",
    )

    parser.add_argument(
        "--depth",
        required=False,
        type=int,
        default=10,
        help="depth to consider for label propagation. [default: 10]",
    )

    parser.add_argument(
        "--mg_threshold",
        required=False,
        type=float,
        default=0.5,
        help="length threshold to consider marker genes. [default: 0.5]",
    )

    parser.add_argument(
        "--bin_mg_threshold",
        required=False,
        type=float,
        default=0.33333,
        help="minimum fraction of marker genes that should be present in a bin. [default: 0.33333]",
    )

    parser.add_argument(
        "--min_bin_size",
        required=False,
        type=int,
        default=200000,
        help="minimum size of a bin to output in base pairs. [default: 200000]",
    )

    parser.add_argument(
        "--delimiter",
        required=False,
        type=str,
        default=",",
        help="delimiter for output results. Supports a comma (,), a semicolon (;), a tab ($'\\t'), a space (\" \") and a pipe (|) [default: , (comma)]",
    )

    parser.add_argument(
        "--nthreads",
        required=False,
        type=int,
        default=8,
        help="number of threads to use. [default: 8]",
    )

    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + version
    )

    args = vars(parser.parse_args())

    return args


def validate(args):
    # Validation of inputs
    # ---------------------------------------------------

    # Check assembler name
    assemblers = ["spades", "megahit", "flye"]
    if args["assembler"].lower() not in assemblers:
        print("\nPlease make sure to provide the correct assembler type.")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Check assembly graph file
    if not os.path.isfile(args["graph"]):
        print("\nFailed to open the assembly graph file.")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Check contigs file
    if not os.path.isfile(args["contigs"]):
        print("\nFailed to open the contigs file.")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Check if paths file is provided when the assembler type is SPAdes
    if args["assembler"].lower() == "spades" and args["paths"] is None:
        print("\nPlease make sure to provide the path to the contigs.paths file.")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Check contigs.paths file for SPAdes
    if args["assembler"].lower() == "spades" and not os.path.isfile(args["paths"]):
        print("\nFailed to open the contigs.paths file.")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Check if paths file is provided when the assembler type is Flye
    if args["assembler"].lower() == "flye" and args["paths"] is None:
        print("\nPlease make sure to provide the path to the assembly_info.txt file.")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Check contigs.paths file for Flye
    if args["assembler"].lower() == "flye" and not os.path.isfile(args["paths"]):
        print("\nFailed to open the assembly_info.txt file.")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Skip paths file when the assembler type is MEGAHIT
    if args["assembler"].lower() == "megahit":
        args["paths"] = "None"

    # Check if abundance file is provided
    if args["abundance"] is None:
        print("\nPlease make sure to provide the path to the abundance file.")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Handle for missing trailing forwardslash in output folder path
    if args["output"][-1:] != "/":
        args["output"] = args["output"] + "/"

    # Create output folder if it does not exist
    if not os.path.isdir(args["output"]):
        subprocess.run("mkdir -p " + args["output"], shell=True)

    # Validate prefix
    if args["prefix"] != "":
        if not args["prefix"].endswith("_"):
            args["prefix"] = args["prefix"] + "_"
    else:
        args["prefix"] = ""

    # Validate min_length
    if args["min_length"] <= 0:
        print("\nPlease enter a valid number for min_length")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Validate p_intra
    if args["p_intra"] <= 0 or args["p_intra"] > 1:
        print("\nPlease enter a valid number for p_intra")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Validate p_inter
    if args["p_inter"] <= 0 or args["p_inter"] > 1:
        print("\nPlease enter a valid number for p_inter")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Validate difference of p_intra and p_inter
    if args["p_inter"] <= 0:
        print(
            "\np_inter cannot be larger than p_intra. Please enter valid numbers for p_intra and p_inter"
        )
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Validate mg_threshold
    if args["mg_threshold"] <= 0 or args["mg_threshold"] > 1:
        print("\nPlease enter a valid number for mg_threshold")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Validate bin_mg_threshold
    if args["bin_mg_threshold"] <= 0 or args["bin_mg_threshold"] > 1:
        print("\nPlease enter a valid number for bin_mg_threshold")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Validate min_bin_size
    if args["min_bin_size"] <= 0:
        print("\nPlease enter a valid number for min_bin_size")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Validate depth
    if args["depth"] <= 0:
        print("\nPlease enter a valid number for depth")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Validate d_limit
    if args["d_limit"] <= 0:
        print("\nPlease enter a valid number for d_limit")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Validate delimiter
    delimiters = [",", ";", " ", "\t", "|"]

    if args["delimiter"] not in delimiters:
        print("\nPlease enter a valid delimiter")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    # Validate number of threads
    if args["nthreads"] <= 0:
        print("\nPlease enter a valid number for the number of threads")
        print("Exiting MetaCoAG...\nBye...!\n")
        sys.exit(1)

    return args
