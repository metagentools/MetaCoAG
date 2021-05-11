#!/usr/bin/env python3

import sys
import csv
import time
import argparse
import re
import heapq
import os
import math
import networkx as nx
import logging
import numpy as np
import operator
import gc
import subprocess
import concurrent.futures

from multiprocessing import Pool
from Bio import SeqIO
from igraph import *
from collections import defaultdict
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
MG_LENGTH_THRESHOLD = 0.3
SEED_MG_THRESHOLD = 0.1


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
ap.add_argument("--w_intra", required=False, type=float, default=2,
                help="maximum weight of an edge matching to assign to the same bin. [default: 2]")
ap.add_argument("--w_inter", required=False, type=int, default=80,
                help="minimum weight of an edge matching to create a new bin. [default: 80]")
ap.add_argument("--depth", required=False, type=int, default=10,
                help="depth to consider for label propagation. [default: 10]")
ap.add_argument("--d_limit", required=False, type=int, default=10,
                help="distance limit for contig matching. [default: 10]")
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
w_intra = args["w_intra"]
w_inter = args["w_inter"]
depth = args["depth"]
d_limit = args["d_limit"]
delimiter = args["delimiter"]
nthreads = args["nthreads"]

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
logger.info("Final binning output file: "+output_path)
logger.info(
    "Minimum length of contigs to consider for compositional probability: "+str(min_length))
logger.info("Depth to consider for label propagation: "+str(depth))
logger.info("w_intra: "+str(w_intra))
logger.info("w_inter: "+str(w_inter))
logger.info("d_limit: "+str(d_limit))
logger.info("Number of threads: "+str(nthreads))

logger.info("MetaCoAG started")

start_time = time.time()


# Get links of the assembly graph
# --------------------------------------------------------

try:
    if assembler == "spades":
        paths, segment_contigs, node_count, contigs_map, contig_names = graph_utils.get_segment_paths_spades(
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
    logger.info("Total number of edges in the assembly graph: " +
                str(len(edge_list)))

    assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

except:
    logger.error(
        "Please make sure that the correct path to the assembly graph file is provided.")
    logger.info("Exiting MetaCoAG... Bye...!")


if assembler == "megahit":

    # Map original contig IDs to contig IDS of assembly graph
    # --------------------------------------------------------

    graph_to_contig_map = BidirectionalMap()

    for (n, m), (n2, m2) in zip(graph_contigs.items(), original_contigs.items()):
        if m == m2:
            graph_to_contig_map[n] = n2

    graph_to_contig_map_rev = graph_to_contig_map.inverse


# Get length and coverage of contigs
# --------------------------------------------------------

logger.info("Obtaining lengths and coverage values of contigs")

if assembler == "megahit":
    seqs, coverages, contig_lengths, zero_cov_contigs = feature_utils.get_cov_len_megahit(
        contigs_file, contig_names_rev, graph_to_contig_map_rev, abundance_file)

else:
    seqs, coverages, contig_lengths, zero_cov_contigs = feature_utils.get_cov_len(
        contigs_file, contig_names_rev, abundance_file)


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
        contigs_file, contig_names_rev, graph_to_contig_map_rev, MG_LENGTH_THRESHOLD, contig_lengths, min_length)

else:
    marker_contigs, marker_contig_counts, contig_markers = marker_gene_utils.get_contigs_with_marker_genes(
        contigs_file, contig_names_rev, MG_LENGTH_THRESHOLD, contig_lengths, min_length)


# Get marker gene counts to make bins
# -----------------------------------------------------

logger.info("Determining seed marker genes")

my_gene_counts = marker_gene_utils.get_seed_marker_gene_counts(
    marker_contig_counts, SEED_MG_THRESHOLD)
my_gene_counts.sort(reverse=True)

# Get the marker genes with maximum count and sort them by the total contig length
total_mg_contig_lengths = {}

for item in marker_contig_counts:

    if marker_contig_counts[item] == my_gene_counts[0]:

        total_contig_lengths = 0

        for contig in marker_contigs[item]:
            length = contig_lengths[contig]
            total_contig_lengths += length

        total_mg_contig_lengths[item] = total_contig_lengths

total_mg_contig_lengths_sorted = sorted(
    total_mg_contig_lengths.items(), key=operator.itemgetter(1), reverse=True)

# Get contigs containing each marker gene for each iteration
seed_iter = {}

n = 0
for item in total_mg_contig_lengths_sorted:
    seed_iter[n] = marker_contigs[item[0]]
    n += 1

for i in range(1, len(my_gene_counts)):
    for item in marker_contig_counts:
        if marker_contig_counts[item] == my_gene_counts[i]:
            seed_iter[n] = marker_contigs[item]
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

logger.info("Number of initial bins detected: "+str(len(seed_iter[0])))

logger.debug("Initialised bins:")
logger.debug(bins)


# Assign contigs with marker genes to bins
# -----------------------------------------------------

edge_weights_per_iteration = {}

logger.info("Matching and assigning contigs with marker genes to bins")

for i in range(len(seed_iter)):

    logger.debug("Iteration "+str(i)+": " +
                 str(len(seed_iter[i]))+" contigs with seed marker genes")

    if i > 0:

        B = nx.Graph()

        common = set(binned_contigs_with_markers).intersection(
            set(seed_iter[i]))

        to_bin = list(set(seed_iter[i]) - common)

        logger.debug(str(len(to_bin))+" contigs to bin in the iteration")

        n_bins = len(bins)

        bottom_nodes = []

        for n in range(n_bins):
            contigid = bins[n][0]
            if contigid not in bottom_nodes:
                bottom_nodes.append(contigid)

        top_nodes = []
        edges = []

        binned_count = 0

        if len(to_bin) != 0:

            for contig in to_bin:

                contigid = contig

                if contigid not in zero_cov_contigs:

                    if contigid not in top_nodes:
                        top_nodes.append(contigid)

                    for b in range(n_bins):

                        log_prob_sum = 0
                        n_contigs = len(bins[b])

                        for j in range(n_contigs):

                            tetramer_dist = matching_utils.get_tetramer_distance(normalized_tetramer_profiles[contigid],
                                                                                 normalized_tetramer_profiles[bins[b][j]])
                            prob_comp = matching_utils.get_comp_probability(
                                tetramer_dist)
                            prob_cov = matching_utils.get_cov_probability(
                                coverages[contigid], coverages[bins[b][j]])

                            prob_product = prob_comp * prob_cov

                            log_prob = 0

                            if prob_product != 0.0:
                                log_prob = - \
                                    (math.log(prob_comp, 10) +
                                     math.log(prob_cov, 10))
                            else:
                                log_prob = MAX_WEIGHT

                            log_prob_sum += log_prob

                        if log_prob_sum != float("inf"):
                            edges.append(
                                (bins[b][0], contigid, log_prob_sum/n_contigs))
                        else:
                            edges.append((bins[b][0], contigid, MAX_WEIGHT))

            B.add_nodes_from(top_nodes, bipartite=0)
            B.add_nodes_from(bottom_nodes, bipartite=1)

            edge_weights = {}

            # Add edges only between nodes of opposite node sets
            for edge in edges:
                edge_weights[(edge[0], edge[1])] = edge[2]
                B.add_edge(edge[0], edge[1], weight=edge[2])

            edge_weights_per_iteration[i] = edge_weights

            top_nodes = {n for n, d in B.nodes(
                data=True) if d['bipartite'] == 0}
            bottom_nodes = set(B) - top_nodes

            if len(top_nodes) > 0:

                my_matching = nx.algorithms.bipartite.matching.minimum_weight_full_matching(
                    B, top_nodes, "weight")

                not_binned = {}

                for l in my_matching:

                    if l in bin_of_contig:

                        b = bin_of_contig[l]

                        if my_matching[l] not in bins[b] and (l, my_matching[l]) in edge_weights:

                            path_len_sum = 0

                            for contig_in_bin in bins[b]:
                                shortest_paths = assembly_graph.get_shortest_paths(
                                    my_matching[l], to=contig_in_bin)

                                if len(shortest_paths) != 0:
                                    path_len_sum += len(shortest_paths[0])

                            avg_path_len = path_len_sum/len(bins[b])

                            logger.debug("Contig with seed MG: " + str(l) + ", Contig to assign: " + str(
                                my_matching[l]) + ", Weight: " + str(edge_weights[(l, my_matching[l])]))

                            if edge_weights[(l, my_matching[l])] <= w_intra and math.floor(avg_path_len) <= d_limit:

                                if len(set(bin_markers[b]).intersection(set(contig_markers[my_matching[l]]))) == 0:

                                    bins[b].append(my_matching[l])
                                    bin_of_contig[my_matching[l]] = b
                                    binned_contigs_with_markers.append(
                                        my_matching[l])
                                    binned_count += 1

                                    bin_markers[b] = list(
                                        set(bin_markers[b] + contig_markers[my_matching[l]]))

                                else:
                                    not_binned[my_matching[l]] = (l, b)

                            else:
                                not_binned[my_matching[l]] = (l, b)

                for nb in not_binned:

                    if edge_weights_per_iteration[i][(not_binned[nb][0], nb)] >= w_inter:

                        path_len_sum = 0

                        for contig_in_bin in bins[not_binned[nb][1]]:

                            shortest_paths = assembly_graph.get_shortest_paths(
                                nb, to=contig_in_bin)

                            if len(shortest_paths) != 0:
                                path_len_sum += len(shortest_paths[0])

                        avg_path_len = path_len_sum / \
                            len(bins[not_binned[nb][1]])

                        if math.floor(avg_path_len) >= d_limit:
                            logger.debug("Creating new bin...")
                            bins[n_bins] = [nb]
                            bin_of_contig[nb] = n_bins
                            binned_count += 1

                            bin_markers[n_bins] = contig_markers[nb]
                            n_bins += 1
                            binned_contigs_with_markers.append(nb)

        logger.debug(str(binned_count)+" contigs binned in the iteration")

logger.debug("Bins with contigs containing seed marker genes")

for b in bins:
    logger.debug(str(b) + ": "+str(bins[b]))

logger.debug("Number of seed MG contigs binned: "+str(len(bin_of_contig)))

if len(seed_iter) > 0:
    del edge_weights_per_iteration
    del B
    del my_matching
    del not_binned
    del edge_weights
    del common
    del to_bin
    del top_nodes
    del bottom_nodes
    del edges

del seed_iter
del my_gene_counts
del marker_contigs
del marker_contig_counts
del total_mg_contig_lengths
gc.collect()


# Further assign contigs with seed marker genes
# -------------------------------------------------

unbinned_mg_contigs = list(
    set(contig_markers.keys()) - set(binned_contigs_with_markers))

unbinned_mg_contig_lengths = {}

for contig in unbinned_mg_contigs:

    contigid = contig

    if contigid not in zero_cov_contigs:
        unbinned_mg_contig_lengths[contig] = contig_lengths[contigid]

unbinned_mg_contig_lengths_sorted = sorted(
    unbinned_mg_contig_lengths.items(), key=operator.itemgetter(1), reverse=True)

logger.debug("Number of unbinned seed MG contigs: " +
             str(len(unbinned_mg_contigs)))

logger.info("Further assign contigs with seed marker genes")

for contig in unbinned_mg_contig_lengths_sorted:

    if contig[1] >= min_length:

        possible_bins = []

        for b in bin_markers:
            common_mgs = list(set(bin_markers[b]).intersection(
                set(contig_markers[contig[0]])))
            if len(common_mgs) == 0:
                possible_bins.append(b)

        if len(possible_bins) != 0:

            contigid = contig[0]

            bin_weights = []

            # bin_path_lens = []

            for b in possible_bins:

                log_prob_sum = 0
                n_contigs = len(bins[b])

                # path_len_sum = 0

                for j in range(n_contigs):

                    tetramer_dist = matching_utils.get_tetramer_distance(normalized_tetramer_profiles[contigid],
                                                                         normalized_tetramer_profiles[bins[b][j]])
                    prob_comp = matching_utils.get_comp_probability(
                        tetramer_dist)
                    prob_cov = matching_utils.get_cov_probability(
                        coverages[contigid], coverages[bins[b][j]])

                    prob_product = prob_comp * prob_cov

                    log_prob = 0

                    if prob_product != 0.0:
                        log_prob = - (math.log(prob_comp, 10) +
                                      math.log(prob_cov, 10))
                    else:
                        log_prob = MAX_WEIGHT

                    log_prob_sum += log_prob

                    # shortest_paths = assembly_graph.get_shortest_paths(contigid, to=bins[b][j])

                    if len(shortest_paths) != 0:
                        path_len_sum += len(shortest_paths[0])

                if log_prob_sum != float("inf"):
                    bin_weights.append(log_prob_sum/n_contigs)
                else:
                    bin_weights.append(MAX_WEIGHT)

            min_b_index = -1

            min_b_index, min_b_value = min(
                enumerate(bin_weights), key=operator.itemgetter(1))

            if min_b_index != -1 and min_b_value <= w_intra:
                bins[possible_bins[min_b_index]].append(contigid)
                bin_of_contig[contigid] = possible_bins[min_b_index]
                binned_contigs_with_markers.append(contigid)

                bin_markers[possible_bins[min_b_index]] = list(set(
                    bin_markers[possible_bins[min_b_index]] + contig_markers[contigid]))

unbinned_mg_contigs = list(
    set(contig_markers.keys()) - set(binned_contigs_with_markers))
logger.debug("Remaining number of unbinned MG seed contigs: " +
             str(len(unbinned_mg_contigs)))
logger.debug(
    "Number of seed MG contigs binned after further assignment: "+str(len(bin_of_contig)))

del unbinned_mg_contigs
del unbinned_mg_contig_lengths
del unbinned_mg_contig_lengths_sorted
gc.collect()


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

logger.info("Propagating labels to unlabelled vertices")

contigs_to_bin = set()

for contig in binned_contigs:
    if contig in non_isolated and contig_lengths[contig] >= min_length:
        closest_neighbours = filter(
            lambda x: x not in binned_contigs and x not in zero_cov_contigs and contig_lengths[x] >= min_length, assembly_graph.neighbors(contig, mode=ALL))
        contigs_to_bin.update(closest_neighbours)

sorted_node_list = []
sorted_node_list_ = [list(label_prop_utils.runBFSLong(x, 1, min_length, binned_contigs, bin_of_contig,
                                                      assembly_graph, normalized_tetramer_profiles, coverages, contig_lengths)) for x in contigs_to_bin]
sorted_node_list_ = [item for sublist in sorted_node_list_ for item in sublist]

for data in sorted_node_list_:
    heapObj = DataWrap(data)
    heapq.heappush(sorted_node_list, heapObj)

while sorted_node_list:
    best_choice = heapq.heappop(sorted_node_list)
    to_bin, binned, bin_, dist, cov_comp_diff = best_choice.data

    has_mg = False

    common_mgs = []

    if to_bin in contig_markers:
        has_mg = True
        common_mgs = list(set(bin_markers[bin_]).intersection(
            set(contig_markers[to_bin])))

    if len(common_mgs) == 0 and to_bin in unbinned_contigs and cov_comp_diff <= w_intra:

        bins[bin_].append(to_bin)
        bin_of_contig[to_bin] = bin_
        binned_contigs.append(to_bin)
        unbinned_contigs.remove(to_bin)

        if has_mg:
            binned_contigs_with_markers.append(to_bin)
            bin_markers[bin_] = list(
                set(bin_markers[bin_] + contig_markers[to_bin]))

        # Discover to_bin's neighbours
        unbinned_neighbours = set(filter(
            lambda x: x not in binned_contigs and x not in zero_cov_contigs and contig_lengths[x] >= min_length, assembly_graph.neighbors(to_bin, mode=ALL)))
        sorted_node_list = list(
            filter(lambda x: x.data[0] not in unbinned_neighbours, sorted_node_list))
        heapq.heapify(sorted_node_list)

        for n in unbinned_neighbours:
            candidates = list(label_prop_utils.runBFSLong(n, 1, min_length, binned_contigs, bin_of_contig,
                                                          assembly_graph, normalized_tetramer_profiles, coverages, contig_lengths))
            for c in candidates:
                heapq.heappush(sorted_node_list, DataWrap(c))


logger.debug("Remaining number of unbinned contigs: " +
             str(len(unbinned_contigs)))
logger.debug("Total number of binned contigs: "+str(len(binned_contigs)))
logger.debug("Number of binned contigs with markers: " +
             str(len(binned_contigs_with_markers)))


# Propagate labels to vertices of unlabelled long contigs
# -----------------------------------------------------

logger.info(
    "Further propagating labels to connected vertices of unlabelled long contigs")

contigs_to_bin = set()

for contig in binned_contigs:
    if contig in non_isolated and contig_lengths[contig] >= min_length:
        closest_neighbours = filter(lambda x: x not in zero_cov_contigs, label_prop_utils.getClosestLongVertices(
            assembly_graph, contig, binned_contigs, contig_lengths, min_length))
        contigs_to_bin.update(closest_neighbours)

sorted_node_list = []
sorted_node_list_ = [list(label_prop_utils.runBFSLong(x, depth, min_length, binned_contigs, bin_of_contig,
                                                      assembly_graph, normalized_tetramer_profiles, coverages, contig_lengths)) for x in contigs_to_bin]
sorted_node_list_ = [item for sublist in sorted_node_list_ for item in sublist]

for data in sorted_node_list_:
    heapObj = DataWrap(data)
    heapq.heappush(sorted_node_list, heapObj)

while sorted_node_list:
    best_choice = heapq.heappop(sorted_node_list)
    to_bin, binned, bin_, dist, cov_comp_diff = best_choice.data

    has_mg = False

    common_mgs = []

    if to_bin in contig_markers:
        has_mg = True
        common_mgs = list(set(bin_markers[bin_]).intersection(
            set(contig_markers[to_bin])))

    if len(common_mgs) == 0 and to_bin in unbinned_contigs and cov_comp_diff <= w_intra:
        bins[bin_].append(to_bin)
        bin_of_contig[to_bin] = bin_
        binned_contigs.append(to_bin)
        unbinned_contigs.remove(to_bin)

        if has_mg:
            binned_contigs_with_markers.append(to_bin)
            bin_markers[bin_] = list(
                set(bin_markers[bin_] + contig_markers[to_bin]))

        # Discover to_bin's neighbours
        unbinned_neighbours = set(filter(lambda x: x not in zero_cov_contigs, label_prop_utils.getClosestLongVertices(
            assembly_graph, to_bin, binned_contigs, contig_lengths, min_length)))
        sorted_node_list = list(
            filter(lambda x: x.data[0] not in unbinned_neighbours, sorted_node_list))
        heapq.heapify(sorted_node_list)

        for n in unbinned_neighbours:
            candidates = list(label_prop_utils.runBFSLong(n, depth, min_length, binned_contigs, bin_of_contig,
                                                          assembly_graph, normalized_tetramer_profiles, coverages, contig_lengths))
            for c in candidates:
                heapq.heappush(sorted_node_list, DataWrap(c))

logger.debug("Remaining number of unbinned contigs: " +
             str(len(unbinned_contigs)))
logger.debug("Total number of binned contigs: "+str(len(binned_contigs)))


# Propagate labels to vertices of unlabelled long contigs in isolated components
# -----------------------------------------------------

logger.info("Further propagating labels to vertices of unlabelled long contigs")

long_unbinned = list(
    filter(lambda contig: contig_lengths[contig] >= min_length, unbinned_contigs))

assigned = [None for itr in long_unbinned]

executor = concurrent.futures.ThreadPoolExecutor(max_workers=nthreads)


def thread_function(n, contig, coverages, normalized_tetramer_profiles, bins, assembly_graph, w_intra, d_limit):
    bin_result = label_prop_utils.assignLong(
        contig, coverages, normalized_tetramer_profiles, bins, assembly_graph, w_intra, d_limit)
    assigned[n] = bin_result


exec_args = []

for n, contig in enumerate(long_unbinned):
    exec_args.append((n, contig, coverages, normalized_tetramer_profiles,
                      bins, assembly_graph, w_intra, d_limit))

for itr in tqdm(executor.map(lambda p: thread_function(*p), exec_args), total=len(long_unbinned)):
    pass

executor.shutdown(wait=True)

put_to_bins = list(filter(lambda x: x is not None, assigned))

if len(put_to_bins) == 0:
    logger.info("No further contigs were binned")

# Assign contigs to bins
for contig, contig_bin in put_to_bins:

    has_mg = False

    common_mgs = []

    if contig in contig_markers:
        has_mg = True
        common_mgs = list(set(bin_markers[contig_bin]).intersection(
            set(contig_markers[contig])))

    if len(common_mgs) == 0:

        bins[contig_bin].append(contig)
        bin_of_contig[contig] = contig_bin
        binned_contigs.append(contig)
        unbinned_contigs.remove(contig)

        if has_mg:
            binned_contigs_with_markers.append(contig)
            bin_markers[contig_bin] = list(
                set(bin_markers[contig_bin] + contig_markers[contig]))

logger.info("Remaining number of unbinned contigs: " +
            str(len(unbinned_contigs)))
logger.info("Total number of binned contigs: "+str(len(binned_contigs)))


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
