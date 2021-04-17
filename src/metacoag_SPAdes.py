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

from multiprocessing import Pool
from Bio import SeqIO
from igraph import *
from collections import defaultdict
from tqdm import tqdm
from functools import partial

from metacoag_utils import feature_utils
from metacoag_utils import marker_gene_utils
from metacoag_utils import matching_utils
from metacoag_utils import label_prop_utils
from metacoag_utils import graph_utils
from metacoag_utils.label_prop_utils import DataWrap
from metacoag_utils.bidirectionalmap import BidirectionalMap

# Set paramters
#---------------------------------------------------

MAX_WEIGHT = sys.float_info.max
MG_LENGTH_THRESHOLD = 0.5
SEED_MG_THRESHOLD = 0.333333
DEPTH = 1


# Setup argument parser
#---------------------------------------------------

ap = argparse.ArgumentParser(description="""MetaCoAG is a NGS data-based metagenomic contig binning tool that makes use of the 
connectivity information found in assembly graphs, apart from the composition and coverage information. 
MetaCoAG makes use of single-copy marker genes along with a graph matching technique and a label propagation technique to bin contigs.""")

ap.add_argument("--contigs", required=True, help="path to the contigs file")
ap.add_argument("--graph", required=True, help="path to the assembly graph file")
ap.add_argument("--paths", required=True, help="path to the contigs.paths file")
ap.add_argument("--abundance", required=True, help="path to the abundance file")
ap.add_argument("--output", required=True, help="path to the output folder")
ap.add_argument("--prefix", required=False, default='', help="prefix for the output file")
ap.add_argument("--min_length", required=False, type=int, default=1000, help="minimum length of contigs to consider for compositional probability. [default: 1000]")
ap.add_argument("--w_intra", required=False, type=float, default=2, help="maximum weight of an edge matching to assign to the same bin. [default: 2]")
ap.add_argument("--w_inter", required=False, type=int, default=80, help="minimum weight of an edge matching to create a new bin. [default: 80]")
ap.add_argument("--d_limit", required=False, type=int, default=10, help="distance limit for contig matching. [default: 10]")
ap.add_argument("--delimiter", required=False, type=str, default=",", help="delimiter for output results. [default: , (comma)]")
ap.add_argument("--nthreads", required=False, type=int, default=8, help="number of threads to use. [default: 8]")

args = vars(ap.parse_args())

contigs_file = args["contigs"]
assembly_graph_file = args["graph"]
contig_paths_file = args["paths"]
abundance_file = args["abundance"]
output_path = args["output"]
prefix = args["prefix"]
min_length = args["min_length"]
w_intra = args["w_intra"]
w_inter = args["w_inter"]
d_limit = args["d_limit"]
delimiter = args["delimiter"]
nthreads = args["nthreads"]

n_bins = 0

# Setup logger
#-----------------------
logger = logging.getLogger('MetaCoaAG 0.1')
logger.setLevel(logging.DEBUG)
logging.captureWarnings(True)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
consoleHeader = logging.StreamHandler()
consoleHeader.setFormatter(formatter)
consoleHeader.setLevel(logging.INFO)
logger.addHandler(consoleHeader)

# Setup output path for log file
#---------------------------------------------------

fileHandler = logging.FileHandler(output_path+"/"+prefix+"metacoag.log")
fileHandler.setLevel(logging.DEBUG)
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)

logger.info("Welcome to MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs.")
logger.info("This version of MetaCoAG makes use of the assembly graph produced by SPAdes which is based on the de Bruijn graph approach.")

logger.info("Input arguments:")
logger.info("Contigs file: "+contigs_file)
logger.info("Assembly graph file: "+assembly_graph_file)
logger.info("Contig paths file: "+contig_paths_file)
logger.info("Final binning output file: "+output_path)
logger.info("Minimum length of contigs to consider for compositional probability: "+str(min_length))
logger.info("Depth: "+str(DEPTH))
logger.info("w_intra: "+str(w_intra))
logger.info("w_inter: "+str(w_inter))
logger.info("d_limit: "+str(d_limit))
logger.info("Number of threads: "+str(nthreads))

logger.info("MetaCoAG started")

start_time = time.time()


# Get segment paths of the assembly graph
#--------------------------------------------------------

try:
    paths, segment_contigs, node_count, my_map, contig_names = graph_utils.get_segment_paths_spades(contig_paths_file)

except:
    logger.error("Please make sure that the correct path to the contig paths file is provided.")
    logger.info("Exiting MetaCoAG... Bye...!")

contigs_map = my_map
contigs_map_rev = my_map.inverse

## Construct the assembly graph
#-------------------------------

try:
    # Create graph
    assembly_graph = Graph()

    # Add vertices
    assembly_graph.add_vertices(node_count)
    logger.info("Total number of contigs available: "+str(node_count))

    # Name vertices
    for i in range(node_count):
        assembly_graph.vs[i]["id"] = i
        assembly_graph.vs[i]["label"] = "NODE_"+str(contigs_map[i])

    # Get list of edges
    edge_list = graph_utils.get_graph_edges_spades(assembly_graph_file, node_count, 
                                                   contigs_map, contigs_map_rev, paths, segment_contigs)

    # Add edges to the graph
    assembly_graph.add_edges(edge_list)
    logger.info("Total number of edges in the assembly graph: "+str(len(edge_list)))

    assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

except:
    logger.error("Please make sure that the correct path to the assembly graph file is provided.")
    logger.info("Exiting MetaCoAG... Bye...!")


# Get length and coverage of contigs
#--------------------------------------------------------

logger.info("Obtaining lengths and coverage values of contigs")

seqs, coverages, contig_lengths = feature_utils.get_cov_len_spades(contigs_file, contigs_map_rev, abundance_file)


# Get tetramer composition of contigs
#--------------------------------------------------------

logger.info("Obtaining tetranucleotide frequencies of contigs")

tetramer_profiles, normalized_tetramer_profiles = feature_utils.get_tetramer_profiles(output_path, seqs, nthreads)


# Get contigs with marker genes
#-----------------------------------------------------

logger.info("Scanning for single-copy marker genes")

if not os.path.exists(contigs_file+".hmmout"):
    logger.info("Obtaining hmmout file")
    marker_gene_utils.scan_for_marker_genes(contigs_file, nthreads)
else:
    logger.info(".hmmout file already exists")


logger.info("Obtaining contigs with single-copy marker genes")

marker_contigs, marker_contig_counts, contig_markers = marker_gene_utils.get_contigs_with_marker_genes(contigs_file, MG_LENGTH_THRESHOLD)


# Get marker gene counts to make bins
#-----------------------------------------------------

logger.info("Determining seed marker genes")

my_gene_counts = marker_gene_utils.get_seed_marker_gene_counts(marker_contig_counts, SEED_MG_THRESHOLD)
my_gene_counts.sort(reverse=True)

# Get the marker genes with maximum count and sort them by the total contig length
total_mg_contig_lengths = {}

for item in marker_contig_counts:

    if marker_contig_counts[item] == my_gene_counts[0]:

        total_contig_lengths = 0

        for contig in marker_contigs[item]:
            start = 'NODE_'
            end = ''
            contig_num = int(re.search('%s(.*)%s' % (start, end), contig).group(1))

            length = contig_lengths[contig_num]
            total_contig_lengths += length

        total_mg_contig_lengths[item] = total_contig_lengths

total_mg_contig_lengths_sorted = sorted(total_mg_contig_lengths.items(),key=operator.itemgetter(1),reverse=True)

# Get contigs containing each marker gene for each iteration
seed_iter = {}

n = 0
for item in total_mg_contig_lengths_sorted:
    seed_iter[n] = marker_contigs[item[0]]
    n += 1

for i in range(1,len(my_gene_counts)):
    for item in marker_contig_counts:
        if marker_contig_counts[item] == my_gene_counts[i]:
            seed_iter[n] = marker_contigs[item]
            n += 1


# Initialise bins
#-----------------------------------------------------

bins = {}
bin_of_contig = {}
bin_markers = {}

binned_contigs_with_markers = []

logger.info("Initialising bins")

for i in range(len(seed_iter[0])):
    
    start = 'NODE_'
    end = ''
    binned_contigs_with_markers.append(seed_iter[0][i])
    contig_num = int(re.search('%s(.*)%s' % (start, end), seed_iter[0][i]).group(1))
    
    bins[i] = [contigs_map_rev[contig_num]]
    bin_of_contig[contigs_map_rev[contig_num]] = i

    bin_markers[i] = contig_markers[seed_iter[0][i]]

logger.info("Number of initial bins detected: "+str(len(seed_iter[0])))

logger.debug("Initialised bins:")
logger.debug(bins)


# Assign contigs with marker genes to bins
#-----------------------------------------------------

edge_weights_per_iteration = {}

logger.info("Matching and assigning contigs with marker genes to bins")
        
for i in range(len(seed_iter)):
    
    logger.debug("Iteration "+str(i)+": "+str(len(seed_iter[i]))+" contigs")
    # print("\nIteration "+str(i)+": "+str(len(seed_iter[i]))+" contigs")
    
    if i>0:
        
        B = nx.Graph()
                
        common = set(binned_contigs_with_markers).intersection(set(seed_iter[i]))
        
        to_bin = list(set(seed_iter[i]) - common)
        
        n_bins = len(bins)
        
        bottom_nodes = []

        for n in range(n_bins):
            contigid = bins[n][0]
            if contigid not in bottom_nodes:
                bottom_nodes.append(contigid)
               
        top_nodes = []
        edges = []
        
        if len(to_bin) != 0:
        
            for contig in to_bin:

                start = 'NODE_'
                end = ''
                contig_num = int(re.search('%s(.*)%s' % (start, end), contig).group(1))

                contigid = contigs_map_rev[contig_num]

                if contigid not in top_nodes:
                    top_nodes.append(contigid)

                for b in range(n_bins):

                    log_prob_sum = 0
                    n_contigs = len(bins[b])

                    for j in range(n_contigs):

                        tetramer_dist = matching_utils.get_tetramer_distance(normalized_tetramer_profiles[contigid], 
                                                                             normalized_tetramer_profiles[bins[b][j]])
                        prob_comp = matching_utils.get_comp_probability(tetramer_dist)
                        prob_cov = matching_utils.get_cov_probability(coverages[contigid], coverages[bins[b][j]])

                        prob_product = prob_comp * prob_cov

                        log_prob = 0

                        if prob_product != 0.0:
                            log_prob = - (math.log(prob_comp, 10) + math.log(prob_cov, 10))
                        else:
                            log_prob = MAX_WEIGHT

                        log_prob_sum += log_prob

                    if log_prob_sum != float("inf"):
                        edges.append((bins[b][0], contigid, log_prob_sum/n_contigs))
                    else:
                        edges.append((bins[b][0], contigid, MAX_WEIGHT))

            B.add_nodes_from(top_nodes, bipartite=0)
            B.add_nodes_from(bottom_nodes, bipartite=1)

            edge_weights = {}


            # Add edges only between nodes of opposite node sets
            for edge in edges:
                edge_weights[(edge[0], edge[1])] = edge[2]
                B.add_edge(edge[0], edge[1], weight = edge[2])

            edge_weights_per_iteration[i] = edge_weights

            top_nodes = {n for n, d in B.nodes(data=True) if d['bipartite'] == 0}
            bottom_nodes = set(B) - top_nodes

            my_matching = nx.algorithms.bipartite.matching.minimum_weight_full_matching(B, top_nodes, "weight")

            not_binned = {}

            for l in my_matching:

                if l in bin_of_contig:

                    b = bin_of_contig[l]

                    if my_matching[l] not in bins[b] and (l, my_matching[l]) in edge_weights:
                        
                        path_len_sum = 0
                    
                        for contig_in_bin in bins[b]:
                            shortest_paths = assembly_graph.get_shortest_paths(my_matching[l], to=contig_in_bin)
                        
                            if len(shortest_paths) != 0:
                                path_len_sum += len(shortest_paths[0])

                        avg_path_len = path_len_sum/len(bins[b])

                        # print(l, my_matching[l], edge_weights[(l, my_matching[l])], avg_path_len)
                                               
                        if edge_weights[(l, my_matching[l])] <= w_intra and math.floor(avg_path_len) <= d_limit:
                            
                            if len(set(bin_markers[b]).intersection(set(contig_markers["NODE_"+str(contigs_map[my_matching[l]])]))) == 0:

                                bins[b].append(my_matching[l])
                                bin_of_contig[my_matching[l]] = b
                                binned_contigs_with_markers.append("NODE_"+str(contigs_map[my_matching[l]]))
                                
                                bin_markers[b] = list(set(bin_markers[b] + contig_markers["NODE_"+str(contigs_map[my_matching[l]])]))

                            else:
                                not_binned[my_matching[l]] = (l, b)

                        else:
                            not_binned[my_matching[l]] = (l, b)

            for nb in not_binned:
                    
                if edge_weights_per_iteration[i][(not_binned[nb][0], nb)] >= w_inter:
                    
                    path_len_sum = 0
                    
                    for contig_in_bin in bins[not_binned[nb][1]]:

                        shortest_paths = assembly_graph.get_shortest_paths(nb, to=contig_in_bin)
                        
                        if len(shortest_paths) != 0:
                            path_len_sum += len(shortest_paths[0])
                    
                    avg_path_len = path_len_sum/len(bins[not_binned[nb][1]])
                    
                    if math.floor(avg_path_len) >= d_limit:
                        logger.debug("Creating new bin...")
                        bins[n_bins] = [nb]
                        bin_of_contig[nb] = n_bins
                        
                        bin_markers[n_bins] = contig_markers["NODE_"+str(contigs_map[nb])]
                        n_bins += 1
                        binned_contigs_with_markers.append("NODE_"+str(contigs_map[nb]))

logger.debug("Bins with contigs containing seed marker genes")

for b in bins:
    logger.debug(str(b)+ ": "+str(bins[b]))


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


# Write intermediate result to output file
#------------------------------------------

output_bins = []

for contig in bin_of_contig:
    line = []
    line.append("NODE_"+str(contigs_map[contig]))
    line.append(bin_of_contig[contig]+1)
    output_bins.append(line)

output_file = output_path + prefix + 'metacoag_seedmg_output.csv'

with open(output_file, mode='w') as output_file:
    output_writer = csv.writer(output_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    for row in output_bins:
        output_writer.writerow(row)


# Further assign contigs with seed marker genes
#-------------------------------------------------

unbinned_mg_contigs = list(set(contig_markers.keys()) - set(binned_contigs_with_markers))

unbinned_mg_contig_lengths = {}

for contig in unbinned_mg_contigs:
    start = 'NODE_'
    end = ''
    contig_num = int(re.search('%s(.*)%s' % (start, end), contig).group(1))

    unbinned_mg_contig_lengths[contig] = contig_lengths[contig_num]

unbinned_mg_contig_lengths_sorted = sorted(unbinned_mg_contig_lengths.items(),key=operator.itemgetter(1),reverse=True)


for contig in unbinned_mg_contig_lengths_sorted:

    if contig[1] >= min_length:

        possible_bins = []

        for b in bin_markers:
            common_mgs = list(set(bin_markers[b]).intersection(set(contig_markers[contig[0]])))
            if len(common_mgs) == 0:
                possible_bins.append(b)

        if len(possible_bins) != 0:

            start = 'NODE_'
            end = ''
            contig_num = int(re.search('%s(.*)%s' % (start, end), contig[0]).group(1))

            contigid = contigs_map_rev[contig_num]

            bin_weights = []

            bin_path_lens = []

            for b in possible_bins:
        
                log_prob_sum = 0
                n_contigs = len(bins[b])

                path_len_sum = 0

                for j in range(n_contigs):

                    tetramer_dist = matching_utils.get_tetramer_distance(normalized_tetramer_profiles[contigid], 
                                                                            normalized_tetramer_profiles[bins[b][j]])
                    prob_comp = matching_utils.get_comp_probability(tetramer_dist)
                    prob_cov = matching_utils.get_cov_probability(coverages[contigid], coverages[bins[b][j]])

                    prob_product = prob_comp * prob_cov

                    log_prob = 0

                    if prob_product != 0.0:
                        log_prob = - (math.log(prob_comp, 10) + math.log(prob_cov, 10))
                    else:
                        log_prob = MAX_WEIGHT

                    log_prob_sum += log_prob

                    shortest_paths = assembly_graph.get_shortest_paths(contigid, to=bins[b][j])
                                
                    if len(shortest_paths) != 0:
                        path_len_sum += len(shortest_paths[0])

                if path_len_sum != 0:
                    bin_path_lens.append(path_len_sum/len(bins[b]))
                else:
                    bin_path_lens.append(MAX_WEIGHT)

                if log_prob_sum != float("inf"):
                    bin_weights.append(log_prob_sum/n_contigs)
                else:
                    bin_weights.append(MAX_WEIGHT)

            min_b_index = -1
            min_dist_index = -1

            min_b_index, min_b_value = min(enumerate(bin_weights), key=operator.itemgetter(1))
            min_dist_index, min_dist_value = min(enumerate(bin_path_lens), key=operator.itemgetter(1))

            if len(assembly_graph.neighbors(contigid)) != 0:
                    
                if min_b_index != -1 and min_dist_index != -1 and min_b_value <= w_intra:
                    bins[possible_bins[min_b_index]].append(contigid)
                    bin_of_contig[contigid] = possible_bins[min_b_index]

                    bin_markers[possible_bins[min_b_index]] = list(set(bin_markers[possible_bins[min_b_index]] + contig_markers["NODE_"+str(contigs_map[contigid])]))

            else:
                if min_b_index != -1 and min_b_value <= w_intra:
                    bins[possible_bins[min_b_index]].append(contigid)
                    bin_of_contig[contigid] = possible_bins[min_b_index]

                    bin_markers[possible_bins[min_b_index]] = list(set(bin_markers[possible_bins[min_b_index]] + contig_markers["NODE_"+str(contigs_map[contigid])]))
        

del contig_markers
del bin_markers
del binned_contigs_with_markers
del unbinned_mg_contigs
del unbinned_mg_contig_lengths
del unbinned_mg_contig_lengths_sorted
gc.collect()


# Get binned and unbinned contigs
#-----------------------------------------------------

binned_contigs = list(bin_of_contig.keys())
    
unbinned_contigs = list(set([x for x in range(node_count)]) - set(binned_contigs))

logger.info("Number of binned contigs: "+str(len(binned_contigs)))
logger.info("Number of unbinned contigs: "+str(len(unbinned_contigs)))


# Get isolated vertices and components without labels
#-----------------------------------------------------

isolated = graph_utils.get_isolated(node_count, assembly_graph)

non_isolated = graph_utils.get_non_isolated(node_count, assembly_graph, binned_contigs)

logger.info("Number of non-isolated contigs: "+str(len(non_isolated)))

non_isolated_unbinned = list(set(non_isolated).intersection(set(unbinned_contigs)))

logger.info("Number of non-isolated unbinned contigs: "+str(len(non_isolated_unbinned)))


# Propagate labels to unlabelled vertices
#-----------------------------------------------------

logger.info("Propagating labels to unlabelled vertices")

# Initialise progress bar
pbar = tqdm(total=len(non_isolated_unbinned))
   
contigs_to_bin = set()

for contig in binned_contigs:
    if contig in non_isolated:
        closest_neighbours = filter(lambda x: x not in binned_contigs, assembly_graph.neighbors(contig, mode=ALL))
        contigs_to_bin.update(closest_neighbours)

sorted_node_list = []
sorted_node_list_ = [list(label_prop_utils.runBFS(x, DEPTH, min_length, binned_contigs, bin_of_contig, 
                                                  assembly_graph, normalized_tetramer_profiles, coverages, contig_lengths)) for x in contigs_to_bin]
sorted_node_list_ = [item for sublist in sorted_node_list_ for item in sublist]

for data in sorted_node_list_:
    heapObj = DataWrap(data)
    heapq.heappush(sorted_node_list, heapObj)

while sorted_node_list:
    best_choice = heapq.heappop(sorted_node_list)    
    to_bin, binned, bin_, dist, cov_comp_diff = best_choice.data
    
    if to_bin in unbinned_contigs:
        bins[bin_].append(to_bin)
        bin_of_contig[to_bin] = bin_
        binned_contigs.append(to_bin)
        unbinned_contigs.remove(to_bin)

        # Update progress bar
        pbar.update(1)
        
        # Discover to_bin's neighbours
        unbinned_neighbours = set(filter(lambda x: x not in binned_contigs, assembly_graph.neighbors(to_bin, mode=ALL)))
        sorted_node_list = list(filter(lambda x: x.data[0] not in unbinned_neighbours, sorted_node_list))
        heapq.heapify(sorted_node_list)
    
        for n in unbinned_neighbours:
            candidates = list(label_prop_utils.runBFS(n, DEPTH, min_length, binned_contigs, bin_of_contig, 
                                                      assembly_graph, normalized_tetramer_profiles, coverages, contig_lengths))
            for c in candidates:
                heapq.heappush(sorted_node_list, DataWrap(c))

# Close progress bar
pbar.close()

logger.info("Remaining number of unbinned contigs: "+str(len(unbinned_contigs)))
logger.info("Total number of binned contigs: "+str(len(binned_contigs)))


# Write intermediate result to output file
#-----------------------------------

output_bins = []

for contig in bin_of_contig:
    line = []
    line.append("NODE_"+str(contigs_map[contig]))
    line.append(bin_of_contig[contig]+1)
    output_bins.append(line)

output_file = output_path + prefix + 'metacoag_intermediate_output.csv'

with open(output_file, mode='w') as output_file:
    output_writer = csv.writer(output_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    for row in output_bins:
        output_writer.writerow(row)


# Bin remaining unbinned contigs
#-----------------------------------

logger.info("Propagating labels to remaining unlabelled vertices")

# Obtain tetramer profiles of contigs in bins
bin_tetramer_profiles = {}

for b in range(len(bins)):

    bin_tetramer_profile = np.zeros(len(tetramer_profiles[0]))

    for j in range(len(bins[b])):

        if contig_lengths[bins[b][j]] >= min_length:
            bin_tetramer_profile = np.add(bin_tetramer_profile, np.array(tetramer_profiles[bins[b][j]]))

    bin_tetramer_profiles[b] = bin_tetramer_profile/max(1, sum(bin_tetramer_profile))

# Assign contigs
with Pool(nthreads) as p:

    assigned = p.map(partial(label_prop_utils.assign, min_length=min_length, contig_lengths=contig_lengths, 
                             normalized_tetramer_profiles=normalized_tetramer_profiles, 
                             bin_tetramer_profiles=bin_tetramer_profiles, n_bins=len(bins)), unbinned_contigs)

p.close()
p.join()

del bin_tetramer_profiles
gc.collect()

put_to_bins = list(filter(lambda x: x is not None, assigned))

if len(put_to_bins) ==  0:
    logger.info("No further contigs were binned")
else:
    logger.info(str(len(put_to_bins))+" contigs were binned")

# Assign contigs to bins
for contig, contig_bin in put_to_bins:
    bins[contig_bin].append(contig)
    bin_of_contig[contig] = contig_bin
    binned_contigs.append(contig)
    unbinned_contigs.remove(contig)

logger.info("Remaining number of unbinned contigs: "+str(len(unbinned_contigs)))
logger.info("Total number of binned contigs: "+str(len(binned_contigs)))


non_isolated = graph_utils.get_non_isolated(node_count, assembly_graph, binned_contigs)

contigs_to_bin = set()

for contig in binned_contigs:
    if contig in non_isolated:
        closest_neighbours = filter(lambda x: x not in binned_contigs, assembly_graph.neighbors(contig, mode=ALL))
        contigs_to_bin.update(closest_neighbours)

sorted_node_list = []
sorted_node_list_ = [list(label_prop_utils.runBFS(x, DEPTH, min_length, binned_contigs, bin_of_contig, 
                                                  assembly_graph, normalized_tetramer_profiles, coverages, contig_lengths)) for x in contigs_to_bin]
sorted_node_list_ = [item for sublist in sorted_node_list_ for item in sublist]

for data in sorted_node_list_:
    heapObj = DataWrap(data)
    heapq.heappush(sorted_node_list, heapObj)

while sorted_node_list:
    best_choice = heapq.heappop(sorted_node_list)    
    to_bin, binned, bin_, dist, cov_comp_diff = best_choice.data
    
    if to_bin in unbinned_contigs:
        bins[bin_].append(to_bin)
        bin_of_contig[to_bin] = bin_
        binned_contigs.append(to_bin)
        unbinned_contigs.remove(to_bin)
        
        # Discover to_bin's neighbours
        unbinned_neighbours = set(filter(lambda x: x not in binned_contigs, assembly_graph.neighbors(to_bin, mode=ALL)))
        sorted_node_list = list(filter(lambda x: x.data[0] not in unbinned_neighbours, sorted_node_list))
        heapq.heapify(sorted_node_list)
    
        for n in unbinned_neighbours:
            candidates = list(label_prop_utils.runBFS(n, DEPTH, min_length, binned_contigs, bin_of_contig, 
                                                      assembly_graph, normalized_tetramer_profiles, coverages, contig_lengths))
            for c in candidates:
                heapq.heappush(sorted_node_list, DataWrap(c))

logger.info("Remaining number of unbinned contigs: "+str(len(unbinned_contigs)))
logger.info("Total number of binned contigs: "+str(len(binned_contigs)))

# Get elapsed time
#-----------------------------------

# Determine elapsed time
elapsed_time = time.time() - start_time

# Print elapsed time for the process
logger.info("Elapsed time: "+str(elapsed_time)+" seconds")


# Write result to output file
#-----------------------------------

output_bins = []

for contig in bin_of_contig:
    line = []
    line.append("NODE_"+str(contigs_map[contig]))
    line.append(bin_of_contig[contig]+1)
    output_bins.append(line)

output_file = output_path + prefix + 'metacoag_output.csv'

with open(output_file, mode='w') as output_file:
    output_writer = csv.writer(output_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    for row in output_bins:
        output_writer.writerow(row)

logger.info("Final binning results can be found at "+str(output_file.name))


# Exit program
#-----------------------------------

logger.info("Thank you for using MetaCoAG!")