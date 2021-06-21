#!/usr/bin/env python3

import sys
import time
import argparse
import os
import math
import logging
import operator
import gc
import subprocess
import concurrent.futures

from Bio import SeqIO
from igraph import *
from tqdm import tqdm

from metacoag_utils import feature_utils
from metacoag_utils import marker_gene_utils
from metacoag_utils import matching_utils
from metacoag_utils import label_prop_utils
from metacoag_utils import graph_utils
from metacoag_utils.label_prop_utils import DataWrap
from metacoag_utils.bidirectionalmap import BidirectionalMap

# Set paramters
# ---------------------------------------------------

MAX_WEIGHT = sys.float_info.max


# Setup argument parser
# ---------------------------------------------------

ap = argparse.ArgumentParser(description="""MetaCoAG is a NGS data-based metagenomic contig binning tool that makes use of the 
connectivity information found in assembly graphs, apart from the composition and coverage information. 
MetaCoAG makes use of single-copy marker genes along with a graph matching technique and a label propagation technique to bin contigs.""")

ap.add_argument("--assembler", required=True,
                help="name of the assembler used")
ap.add_argument("--contigs", required=True, help="path to the contigs file")
ap.add_argument("--graph", required=True,
                help="path to the assembly graph file")
ap.add_argument("--paths", required=True,
                help="path to the contigs.paths file")
ap.add_argument("--abundance", required=True,
                help="path to the abundance file")
ap.add_argument("--output", required=True, help="path to the output folder")
ap.add_argument("--prefix", required=False, default='',
                help="prefix for the output file")
ap.add_argument("--min_length", required=False, type=int, default=1000,
                help="minimum length of contigs to consider for compositional probability. [default: 1000]")
ap.add_argument("--p_intra", required=False, type=float, default=0.1,
                help="minimum probability of an edge matching to assign to the same bin. [default: 0.1]")
ap.add_argument("--p_inter", required=False, type=float, default=0.01,
                help="maximum probability of an edge matching to create a new bin. [default: 0.01]")
ap.add_argument("--depth", required=False, type=int, default=10,
                help="depth to consider for label propagation. [default: 10]")
ap.add_argument("--mg_threshold", required=False, type=float, default=0.5,
                help="length threshold to consider marker genes. [default: 0.5]")
ap.add_argument("--seed_threshold", required=False, type=float, default=0.1,
                help="threshold to consider contigs with seed marker genes. [default: 0.1]")
ap.add_argument("--d_limit", required=False, type=int, default=20,
                help="distance limit for contig matching. [default: 20]")
ap.add_argument("--break_step", required=False, type=float, default=0.0005,
                help="Probability step to stop breaking bins. [default: 0.0005]")
ap.add_argument("--delimiter", required=False, type=str, default=",",
                help="delimiter for output results. [default: , (comma)]")
ap.add_argument("--nthreads", required=False, type=int, default=8,
                help="number of threads to use. [default: 8]")

args = vars(ap.parse_args())

assembler = args["assembler"]
contigs_file = args["contigs"]
assembly_graph_file = args["graph"]
contig_paths_file = args["paths"]
abundance_file = args["abundance"]
output_path = args["output"]
prefix = args["prefix"]
min_length = args["min_length"]
p_intra = args["p_intra"]
p_inter = args["p_inter"]
depth = args["depth"]
mg_threshold = args["mg_threshold"]
seed_threshold = args["seed_threshold"]
d_limit = args["d_limit"]
break_step = args["break_step"]
delimiter = args["delimiter"]
nthreads = args["nthreads"]

bin_threshold = -math.log(p_intra, 10)
break_threshold = -math.log(p_inter, 10)

n_bins = 0

# Setup logger
# -----------------------
logger = logging.getLogger('MetaCoaAG 0.1')
logger.setLevel(logging.DEBUG)
logging.captureWarnings(True)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
consoleHeader = logging.StreamHandler()
consoleHeader.setFormatter(formatter)
consoleHeader.setLevel(logging.INFO)
logger.addHandler(consoleHeader)

# Setup output path for log file
# ---------------------------------------------------

fileHandler = logging.FileHandler(output_path+"/"+prefix+"metacoag.log")
fileHandler.setLevel(logging.DEBUG)
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)

logger.info(
    "Welcome to MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs.")
logger.info("This version of MetaCoAG makes use of the assembly graph produced by SPAdes which is based on the de Bruijn graph approach.")

logger.info("Input arguments:")
logger.info("Assembler used: "+assembler)
logger.info("Contigs file: "+contigs_file)
logger.info("Assembly graph file: "+assembly_graph_file)
logger.info("Contig paths file: "+contig_paths_file)
logger.info("Abundance file: "+abundance_file)
logger.info("Final binning output file: "+output_path)
logger.info(
    "Minimum length of contigs to consider for compositional probability: "+str(min_length))
logger.info("Depth to consider for label propagation: "+str(depth))
logger.info("p_intra: "+str(p_intra))
logger.info("p_inter: "+str(p_inter))
logger.debug("bin_threshold: "+str(bin_threshold))
logger.debug("break_threshold: "+str(break_threshold))
logger.info("mg_threshold: "+str(mg_threshold))
logger.info("seed_threshold: "+str(seed_threshold))
logger.info("d_limit: "+str(d_limit))
logger.info("depth: "+str(depth))
logger.info("Number of threads: "+str(nthreads))

logger.info("MetaCoAG started")

start_time = time.time()


# Get links of the assembly graph
# --------------------------------------------------------

try:
    if assembler == "spades":
        paths, segment_contigs, contig_segments, node_count, contigs_map, contig_names = graph_utils.get_segment_paths_spades(
            contig_paths_file)
        contigs_map_rev = contigs_map.inverse
        contig_names_rev = contig_names.inverse

    if assembler == "megahit":

        original_contigs = {}
        contig_descriptions = {}

        for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
            original_contigs[record.id] = str(record.seq)
            contig_descriptions[record.id] = record.description

        node_count, graph_contigs, links, contig_names = graph_utils.get_links_megahit(
            assembly_graph_file)

        contig_names_rev = contig_names.inverse

    if assembler == "flye":
        node_count, links, contig_names = graph_utils.get_links_flye(
            assembly_graph_file)
        contig_names_rev = contig_names.inverse

except:
    logger.error(
        "Please make sure that the correct path to the contig paths file is provided.")
    logger.info("Exiting MetaCoAG... Bye...!")


# Construct the assembly graph
# -------------------------------

try:
    # Create graph
    assembly_graph = Graph()

    # Add vertices
    assembly_graph.add_vertices(node_count)
    logger.info("Total number of contigs available: "+str(node_count))

    # Name vertices
    for i in range(node_count):
        assembly_graph.vs[i]["id"] = i
        assembly_graph.vs[i]["label"] = contig_names[i]

    # Get list of edges
    if assembler == "spades":
        edge_list = graph_utils.get_graph_edges_spades(assembly_graph_file, node_count,
                                                       contigs_map, contigs_map_rev, paths, segment_contigs)

    if assembler == "flye":
        edge_list = graph_utils.get_graph_edges_flye(links, contig_names_rev)

    if assembler == "megahit":
        edge_list = graph_utils.get_graph_edges_megahit(
            links, contig_names_rev)

    # Add edges to the graph
    assembly_graph.add_edges(edge_list)

    assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

    logger.info("Total number of edges in the assembly graph: " +
                str(len(list(assembly_graph.es))))

except:
    logger.error(
        "Please make sure that the correct path to the assembly graph file is provided.")
    logger.info("Exiting MetaCoAG... Bye...!")


if assembler == "megahit":

    # Map original contig IDs to contig IDS of assembly graph
    graph_to_contig_map = BidirectionalMap()

    for (n, m), (n2, m2) in zip(graph_contigs.items(), original_contigs.items()):
        if m == m2:
            graph_to_contig_map[n] = n2

    graph_to_contig_map_rev = graph_to_contig_map.inverse


# Get length and coverage of contigs
# --------------------------------------------------------

logger.info("Obtaining lengths and coverage values of contigs")

if assembler == "megahit":
    seqs, coverages, contig_lengths, n_samples = feature_utils.get_cov_len_megahit(
        contigs_file, contig_names_rev, graph_to_contig_map_rev, abundance_file, min_length)

else:
    seqs, coverages, contig_lengths, n_samples = feature_utils.get_cov_len(
        contigs_file, contig_names_rev, abundance_file, min_length)


# Set intra weight and inter weight
# --------------------------------------------------------

w_intra = bin_threshold * (n_samples+1)
w_inter = break_threshold * (n_samples+1)

logger.debug("w_intra: "+str(w_intra))
logger.debug("w_inter: "+str(w_inter))


# Get tetramer composition of contigs
# --------------------------------------------------------

logger.info("Obtaining tetranucleotide frequencies of contigs")

tetramer_profiles, normalized_tetramer_profiles = feature_utils.get_tetramer_profiles(
    output_path, seqs, nthreads)


# Get contigs with marker genes
# -----------------------------------------------------

logger.info("Scanning for single-copy marker genes")

if not os.path.exists(contigs_file+".hmmout"):
    logger.info("Obtaining hmmout file")
    marker_gene_utils.scan_for_marker_genes(contigs_file, nthreads)
else:
    logger.info(".hmmout file already exists")


logger.info("Obtaining contigs with single-copy marker genes")

if assembler == "megahit":
    marker_contigs, marker_contig_counts, contig_markers = marker_gene_utils.get_contigs_with_marker_genes_megahit(
        contigs_file, contig_names_rev, graph_to_contig_map_rev, mg_threshold, contig_lengths, min_length)

else:
    marker_contigs, marker_contig_counts, contig_markers = marker_gene_utils.get_contigs_with_marker_genes(
        contigs_file, contig_names_rev, mg_threshold, contig_lengths, min_length)

    all_mg_contigs = marker_gene_utils.get_all_contigs_with_marker_genes(
        contigs_file, contig_names_rev)


# Get marker gene counts to make bins
# -----------------------------------------------------

logger.info("Determining contig counts for each single-copy marker gene")

my_gene_counts = list(marker_contig_counts.values())
my_gene_counts.sort(reverse=True)

logger.debug("Contig counts of single-copy marker genes:")
logger.debug(str(my_gene_counts))


# Get contigs containing each single-copy marker gene for each iteration
# -----------------------------------------------------

seed_iter = {}

n = 0

unique_my_gene_counts = list(set(my_gene_counts))
unique_my_gene_counts.sort(reverse=True)

for g_count in unique_my_gene_counts:

    # Get the marker genes with maximum count and sort them by the total contig length
    total_contig_mgs = {}

    for item in marker_contig_counts:

        if marker_contig_counts[item] == g_count:

            total_contig_lengths = 0

            for contig in marker_contigs[item]:
                contig_mg_counts = len(contig_markers[contig])
                total_contig_lengths += contig_mg_counts

            total_contig_mgs[item] = total_contig_lengths

    total_contig_mgs_sorted = sorted(
        total_contig_mgs.items(), key=operator.itemgetter(1), reverse=True)

    for item in total_contig_mgs_sorted:
        seed_iter[n] = marker_contigs[item[0]]
        n += 1


# Initialise bins
# -----------------------------------------------------

bins = {}
bin_of_contig = {}
bin_markers = {}

binned_contigs_with_markers = []

logger.info("Initialising bins")

for i in range(len(seed_iter[0])):

    binned_contigs_with_markers.append(seed_iter[0][i])
    contig_num = seed_iter[0][i]

    bins[i] = [contig_num]
    bin_of_contig[contig_num] = i

    bin_markers[i] = contig_markers[contig_num]

logger.debug("Number of initial bins detected: "+str(len(seed_iter[0])))

logger.debug("Initialised bins:")
logger.debug(bins)


# Assign contigs with marker genes to bins
# -----------------------------------------------------

logger.info("Matching and assigning contigs with single-copy marker genes to bins")

bins, bin_of_contig, n_bins, bin_markers, binned_contigs_with_markers = matching_utils.match_contigs(seed_iter, bins, n_bins, bin_of_contig, binned_contigs_with_markers, bin_markers, contig_markers, contig_lengths, normalized_tetramer_profiles, coverages, assembly_graph, w_intra, w_inter, d_limit)

logger.debug("Number of bins after matching: " + str(len(bins)))

logger.debug("Bins with contigs containing seed marker genes")

for b in bins:
    logger.debug(str(b) + ": "+str(bins[b]))

logger.debug("Number of binned contigs with single-copy marker genes: "+str(len(bin_of_contig)))

del seed_iter
del my_gene_counts
del marker_contigs
del marker_contig_counts
del total_contig_mgs
gc.collect()


# Further assign contigs with seed marker genes
# -------------------------------------------------

unbinned_mg_contigs = list(
    set(contig_markers.keys()) - set(binned_contigs_with_markers))

unbinned_mg_contig_lengths = {}

for contig in unbinned_mg_contigs:
    contigid = contig
    unbinned_mg_contig_lengths[contig] = contig_lengths[contigid]

unbinned_mg_contig_lengths_sorted = sorted(
    unbinned_mg_contig_lengths.items(), key=operator.itemgetter(1), reverse=True)

logger.debug("Number of unbinned contigs with single-copy marker genes: " +
             str(len(unbinned_mg_contigs)))

logger.info("Further assigning contigs with single-copy marker genes")

bins, bin_of_contig, n_bins, bin_markers, binned_contigs_with_markers = matching_utils.further_match_contigs(unbinned_mg_contig_lengths_sorted, min_length, bins, n_bins, bin_of_contig, binned_contigs_with_markers, bin_markers, contig_markers, normalized_tetramer_profiles, coverages, w_intra)

unbinned_mg_contigs = list(set(contig_markers.keys()) - set(binned_contigs_with_markers))

logger.debug("Remaining number of unbinned MG seed contigs: " + str(len(unbinned_mg_contigs)))
logger.debug("Number of binned contigs with single-copy marker genes: "+str(len(bin_of_contig)))

del unbinned_mg_contigs
del unbinned_mg_contig_lengths
del unbinned_mg_contig_lengths_sorted
gc.collect()


# Get seed bin counts and profiles
# -----------------------------------------------------

seed_bin_count = []

for i in bins:
    seed_bin_count.append(len(bins[i]))

bin_seed_tetramer_profiles, bin_seed_coverage_profiles = feature_utils.get_bin_profiles(bins, coverages, normalized_tetramer_profiles)


# Get binned and unbinned contigs
# -----------------------------------------------------

binned_contigs = list(bin_of_contig.keys())

unbinned_contigs = list(
    set([x for x in range(node_count)]) - set(binned_contigs))

logger.debug("Number of binned contigs: "+str(len(binned_contigs)))
logger.debug("Number of unbinned contigs: "+str(len(unbinned_contigs)))


# Get isolated vertices and components without labels
# -----------------------------------------------------

isolated = graph_utils.get_isolated(node_count, assembly_graph)

non_isolated = graph_utils.get_non_isolated(
    node_count, assembly_graph, binned_contigs)

logger.debug("Number of non-isolated contigs: "+str(len(non_isolated)))

non_isolated_unbinned = list(
    set(non_isolated).intersection(set(unbinned_contigs)))

non_isolated_long_unbinned = list(filter(
    lambda contig: contig_lengths[contig] >= min_length, non_isolated_unbinned))

logger.debug("Number of non-isolated unbinned contigs: " +
            str(len(non_isolated_unbinned)))
logger.debug("Number of non-isolated long unbinned contigs: " +
            str(len(non_isolated_long_unbinned)))
logger.debug("Number of binned contigs with markers: " +
            str(len(binned_contigs_with_markers)))



# Propagate labels to vertices of unlabelled long contigs
# -----------------------------------------------------

logger.info("Propagating labels to connected vertices of unlabelled long contigs")

bins, bin_of_contig, bin_markers, binned_contigs_with_markers = label_prop_utils.label_prop(bin_of_contig, bins, contig_markers, bin_markers, binned_contigs_with_markers, seed_bin_count, non_isolated, contig_lengths, min_length, assembly_graph, normalized_tetramer_profiles, coverages, 1, w_intra)

logger.debug("Total number of binned contigs: "+str(len(bin_of_contig)))


# Further propagate labels to vertices of unlabelled long contigs
# ----------------------------------------------------------------

bins, bin_of_contig, bin_markers, binned_contigs_with_markers = label_prop_utils.label_prop(bin_of_contig, bins, contig_markers, bin_markers, binned_contigs_with_markers, seed_bin_count, non_isolated, contig_lengths, min_length, assembly_graph, normalized_tetramer_profiles, coverages, depth, w_inter)

logger.debug("Total number of binned contigs: "+str(len(bin_of_contig)))


# Get binned and unbinned contigs
# -----------------------------------------------------

binned_contigs = list(bin_of_contig.keys())

unbinned_contigs = list(
    set([x for x in range(node_count)]) - set(binned_contigs))

logger.debug("Number of binned contigs: "+str(len(binned_contigs)))
logger.debug("Number of unbinned contigs: "+str(len(unbinned_contigs)))



# Propagate labels to vertices of unlabelled long contigs in isolated components
# -------------------------------------------------------------------------------

logger.info("Further propagating labels to vertices of unlabelled long contigs")

long_unbinned = list(
    filter(lambda contig: contig not in bin_of_contig and contig_lengths[contig] >= min_length, [x for x in range(node_count)]))

assigned = [None for itr in long_unbinned]

executor = concurrent.futures.ThreadPoolExecutor(max_workers=nthreads)

def thread_function(n, contig, coverages, normalized_tetramer_profiles, bin_seed_tetramer_profiles, bin_seed_coverage_profiles):
    bin_result = label_prop_utils.assign_long(
        contig, coverages, normalized_tetramer_profiles, bin_seed_tetramer_profiles, bin_seed_coverage_profiles)
    assigned[n] = bin_result

exec_args = []

for n, contig in enumerate(long_unbinned):
    exec_args.append(
        (n, contig, coverages, normalized_tetramer_profiles, bin_seed_tetramer_profiles, bin_seed_coverage_profiles))

for itr in tqdm(executor.map(lambda p: thread_function(*p), exec_args), total=len(long_unbinned)):
    pass

executor.shutdown(wait=True)

put_to_bins = [x for x in assigned if x is not None]

if len(put_to_bins) == 0:
    logger.debug("No further contigs were binned")
else:
    bins, bin_of_contig, bin_markers, binned_contigs_with_markers = label_prop_utils.assign_to_bins(put_to_bins, bins, bin_of_contig, bin_markers, binned_contigs_with_markers, contig_markers)

logger.debug("Total number of binned contigs: "+str(len(bin_of_contig)))


#  Further propagate labels to vertices of unlabelled long contigs
# ----------------------------------------------------------------

logger.info("Further propagating labels to connected vertices of unlabelled long contigs")

bins, bin_of_contig, bin_markers, binned_contigs_with_markers = label_prop_utils.final_label_prop(bin_of_contig, bins, contig_markers, bin_markers, binned_contigs_with_markers, seed_bin_count, contig_lengths, min_length, assembly_graph, normalized_tetramer_profiles, coverages, depth, MAX_WEIGHT)

logger.debug("Total number of binned contigs: "+str(len(bin_of_contig)))


# Get elapsed time
# -----------------------------------

# Determine elapsed time
elapsed_time = time.time() - start_time

# Print elapsed time for the process
logger.info("Elapsed time: "+str(elapsed_time)+" seconds")


# Write result to output file
# -----------------------------------

logger.info("Writing the Final Binning result to file")

output_bins_path = output_path + prefix + "bins/"

if not os.path.isdir(output_bins_path):
    subprocess.run("mkdir -p "+output_bins_path, shell=True)

for b in range(len(bins)):

    with open(output_bins_path + "bin_" + str(b+1) + "_ids.txt", "w") as bin_file:
        for contig in bins[b]:

            if assembler == "megahit":
                bin_file.write(
                    contig_descriptions[graph_to_contig_map[contig_names[contig]]]+"\n")
            else:
                bin_file.write(contig_names[contig]+"\n")

    subprocess.run("awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' " + output_bins_path + "bin_" + str(
        b+1) + "_ids.txt " + contigs_file + " > " + output_bins_path + "bin_" + str(b+1) + "_seqs.fasta", shell=True)

logger.info("Final binning results can be found in "+str(output_bins_path))


# Exit program
# -----------------------------------

logger.info("Thank you for using MetaCoAG!")
