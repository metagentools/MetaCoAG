#!/usr/bin/env python3

import sys
import csv
import time
import argparse
import re
import heapq
import os
import itertools
import math
import numpy as np
import networkx as nx
import pathlib

from multiprocessing import Pool
from Bio import SeqIO
from igraph import *
from collections import defaultdict
from scipy.spatial import distance
from scipy.stats import poisson
from tqdm import tqdm
from bidirectionalmap.bidirectionalmap import BidirectionalMap

# Set paramters
#---------------------------------------------------

mu_intra, sigma_intra = 0, 0.01037897 / 2
mu_inter, sigma_inter = 0.0676654, 0.03419337
max_weight = sys.float_info.max
mg_length_threshold = 0.5
seed_mg_threshold = 0.333


# Setup argument parser
#---------------------------------------------------

ap = argparse.ArgumentParser(description="""MetaCoAG is a NGS data-based metagenomic contig binning tool that makes use of the 
connectivity information found in assembly graphs, apart from the composition and coverage information. 
MetaCoAG makes use of single-copy marker genes along with a graph matching technique and a label propagation technique to bin contigs.""")

ap.add_argument("--contigs", required=True, help="path to the contigs file")
ap.add_argument("--graph", required=True, help="path to the assembly graph file")
ap.add_argument("--abundance", required=True, help="path to the abundance file")
ap.add_argument("--output", required=True, help="path to the output folder")
ap.add_argument("--prefix", required=False, default='', help="prefix for the output file")
ap.add_argument("--depth", required=False, type=int, default=5, help="maximum depth for the breadth-first-search. [default: 5]")
ap.add_argument("--min_length", required=False, type=int, default=1000, help="minimum length of contigs to consider for compositional probability. [default: 1000]")
ap.add_argument("--alpha_intra", required=False, type=int, default=2, help="maximum weight of an edge matching to assign to the same bin. [default: 2]")
ap.add_argument("--alpha_inter", required=False, type=int, default=80, help="minimum weight of an edge matching to create a new bin. [default: 80]")
ap.add_argument("--dist_intra", required=False, type=int, default=10, help="maximum distance of a contig matched to assign to the same bin. [default: 10]")
ap.add_argument("--dist_inter", required=False, type=int, default=10, help="minimum distance of a contig matched to create a new bin. [default: 10]")
ap.add_argument("--nthreads", required=False, type=int, default=8, help="number of threads to use. [default: 8]")

args = vars(ap.parse_args())

contigs_file = args["contigs"]
assembly_graph_file = args["graph"]
abundance_file = args["abundance"]
output_path = args["output"]
prefix = args["prefix"]
depth = args["depth"]
min_length = args["min_length"]
alpha_intra = args["alpha_intra"]
alpha_inter = args["alpha_inter"]
dist_intra = args["dist_intra"]
dist_inter = args["dist_inter"]
nthreads = args["nthreads"]

n_bins = 0

print("\nWelcome to MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs.")
print("This version of MetaCoAG makes use of the assembly graph produced by MEGAHIT which is based on the de Bruijn graph approach.")

print("\nInput arguments:", contigs_file)
print("Assembly graph file:", assembly_graph_file)
print("Final binning output file:", output_path)
print("Minimum length of contigs to consider for compositional probability:", min_length)
print("Depth:", depth)
print("alpha_intra:", alpha_intra)
print("alpha_inter:", alpha_inter)
print("dist_intra:", dist_intra)
print("dist_inter:", dist_inter)
print("Number of threads:", nthreads)

print("\nMetaCoAG started\n-------------------")

start_time = time.time()


# Get original contig IDs
#-------------------------------

original_contigs = {}

for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
    original_contigs[record.id] = str(record.seq)


# Build the assembly graph
#--------------------------------------------------------

node_count = 0

graph_contigs = {}

links = []

my_map = BidirectionalMap()

try:

    # Get links from .gfa file
    with open(assembly_graph_file) as file:

        line = file.readline()

        while line != "":

            # Identify lines with link information
            if line.startswith("L"):
                link = []

                strings = line.split("\t")

                start_1 = 'NODE_'
                end_1 = '_length'

                link1 = int(re.search('%s(.*)%s' % (start_1, end_1), strings[1]).group(1))

                start_2 = 'NODE_'
                end_2 = '_length'

                link2 = int(re.search('%s(.*)%s' % (start_2, end_2), strings[3]).group(1))

                link.append(link1)
                link.append(link2)
                links.append(link)

            elif line.startswith("S"):
                strings = line.split()

                start = 'NODE_'
                end = '_length'

                contig_num = int(re.search('%s(.*)%s' % (start, end), strings[1]).group(1))

                my_map[node_count] = int(contig_num)

                graph_contigs[contig_num] = strings[2]

                node_count += 1

            line = file.readline()


    print("Total number of contigs available: "+str(node_count))

    contigs_map = my_map
    contigs_map_rev = my_map.inverse

    # Create graph
    assembly_graph = Graph()

    # Add vertices
    assembly_graph.add_vertices(node_count)

    # Create list of edges
    edge_list = []

    for i in range(node_count):
        assembly_graph.vs[i]["id"]= i
        assembly_graph.vs[i]["label"]= str(contigs_map[i])

    # Iterate links
    for link in links:
        # Remove self loops
        if link[0] != link[1]:
            # Add edge to list of edges
            edge_list.append((contigs_map_rev[link[0]], contigs_map_rev[link[1]]))

    # Add edges to the graph
    assembly_graph.add_edges(edge_list)
    assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

except:
    print("Please make sure that the correct path to the assembly graph file is provided.")
    print("Exiting MetaCoAG... Bye...!")
    sys.exit(1)

print("Total number of edges in the assembly graph: "+str(len(edge_list)))


# Map original contig IDs to contig IDS of assembly graph
#--------------------------------------------------------

graph_to_contig_map = BidirectionalMap()    

for (n,m), (n2,m2) in zip(graph_contigs.items(), original_contigs.items()):
    if m==m2:
        graph_to_contig_map[n] = n2

graph_to_contig_map_rev = graph_to_contig_map.inverse



# Get length and coverage of contigs
#--------------------------------------------------------

print("\nObtaining lengths and coverage values of contigs...")

contig_lengths = {}

longer_than_min_length = 0

i = 0

seqs = []

for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
    
    contig_num = contigs_map_rev[int(graph_to_contig_map_rev[record.id])]
    
    length = len(record.seq)

    if length >= min_length:
        longer_than_min_length += 1
    
    contig_lengths[contig_num] = length
    
    seqs.append(str(record.seq))
    
    i+=1

print("\nNumber of contigs longer than minimum length threshold:", longer_than_min_length)

coverages = {}

with open(abundance_file, "r") as my_file:
    line = my_file.readline().strip()
    
    while line!="":
        strings = line.split("\t")

        contig_num = contigs_map_rev[int(graph_to_contig_map_rev[strings[0]])]

        coverages[contig_num] = int(float(strings[1]))
        
        line = my_file.readline()


# Get tetramer composition of contigs
#--------------------------------------------------------

print("\nObtaining tetranucleotide frequencies of contigs...")

def get_rc(seq):
    rev = reversed(seq)
    return "".join([complements.get(i,i) for i in rev])


def mer2bits(kmer):
    bit_mer=nt_bits[kmer[0]]
    for c in kmer[1:]:
        bit_mer = (bit_mer << 2) | nt_bits[c]
    return bit_mer

def compute_kmer_inds(k):
    kmer_inds = {}
    kmer_count_len = 0

    alphabet = 'ACGT'
    
    all_kmers = [''.join(kmer) for kmer in itertools.product(alphabet,repeat=k)]
    all_kmers.sort()
    ind = 0
    for kmer in all_kmers:
        bit_mer = mer2bits(kmer)
        rc_bit_mer = mer2bits(get_rc(kmer))
        if rc_bit_mer in kmer_inds:
            kmer_inds[bit_mer] = kmer_inds[rc_bit_mer]
        else:
            kmer_inds[bit_mer] = ind
            kmer_count_len += 1
            ind += 1
            
    return kmer_inds, kmer_count_len

def count_kmers(args):
    seq, k, kmer_inds, kmer_count_len = args
    profile = np.zeros(kmer_count_len)
    arrs = []
    seq = list(seq.strip())
    
    for i in range(0, len(seq) - k + 1):
        bit_mer = mer2bits(seq[i:i+k])
        index = kmer_inds[bit_mer]
        profile[index] += 1
    profile = profile/max(1, sum(profile))
    
    return profile


tetramer_profiles = {}

i=0
if os.path.isfile(output_path+"contig_tetramers.txt"):
    with open(output_path+"contig_tetramers.txt") as tetramers_file:
        for line in tetramers_file.readlines():
            f_list = [float(i) for i in line.split(" ") if i.strip()]
            tetramer_profiles[i] = f_list
            i+=1

else:
    complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    nt_bits = {'A':0,'C':1,'G':2,'T':3}

    kmer_inds_4, kmer_count_len_4 = compute_kmer_inds(4)

    pool = Pool(nthreads)
    record_tetramers = pool.map(count_kmers, [(seq, 4, kmer_inds_4, kmer_count_len_4) for seq in seqs])
    pool.close()
    
    i=0

    for l in range(len(record_tetramers)):
        tetramer_profiles[i] = record_tetramers[l]
        i+=1
    
    with open(output_path+"contig_tetramers.txt", "w+") as myfile:
        for l in range(len(record_tetramers)):
            for j in range(len(record_tetramers[l])):
                myfile.write(str(record_tetramers[l][j])+" ")
            myfile.write("\n")


# Get contigs with marker genes
#-----------------------------------------------------

def scan_for_marker_genes(contig_file, hard=0):
# Modified from SolidBin
    
    software_path = pathlib.Path(__file__).parent.absolute()
    
    fragScanURL = os.path.join(software_path.parent, 'auxiliary', 'FragGeneScan1.31', 'run_FragGeneScan.pl')
    hmmExeURL = os.path.join(software_path.parent, 'auxiliary', 'hmmer-3.3', 'src', 'hmmsearch')
    markerURL = os.path.join(software_path.parent, 'auxiliary', 'marker.hmm')

    print(markerURL)
    
    fragResultURL = contig_file+".frag.faa"
    hmmResultURL = contig_file+".hmmout"
    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL+" -genome="+contig_file+" -out="+contig_file + \
            ".frag -complete=0 -train=complete -thread="+str(nthreads)+" 1>" + \
            contig_file+".frag.out 2>"+contig_file+".frag.err"
        print("exec cmd: "+fragCmd)
        os.system(fragCmd)

    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL+" --domtblout "+hmmResultURL+" --cut_tc --cpu "+str(nthreads)+" " + \
                markerURL+" "+fragResultURL+" 1>"+hmmResultURL+".out 2>"+hmmResultURL+".err"
            print("exec cmd: "+hmmCmd)
            os.system(hmmCmd)

        else:
            print("HMMER search failed! Path: "+hmmResultURL + " does not exist.")
    else:
        print("FragGeneScan failed! Path: "+fragResultURL + " does not exist.")

print("\nScanning for single-copy marker genes...")

if not os.path.exists(contigs_file+".hmmout"):
    print("Obtaining hmmout file...")
    scan_for_marker_genes(contigs_file)
else:
    print("hmmout file already exists...")


marker_contigs = {}
marker_contig_counts = {}

print("\nObtaining contigs with single-copy marker genes...")

with open(contigs_file+".hmmout", "r") as myfile:
    for line in myfile.readlines():
        if not line.startswith("#"):
            strings = line.strip().split()

            name_strings = strings[0].split("_")
            
            contig = strings[0]
            marker_gene = strings[3]
            marker_gene_length = int(strings[5])
            
            mapped_marker_length = int(strings[16]) - int(strings[15])
            
            contig_num = '_'.join(name_strings[0:2])
            
            if mapped_marker_length > marker_gene_length*mg_length_threshold:
                
                if marker_gene not in marker_contigs:
                    marker_contigs[marker_gene] = [contig_num]
                else:
                    marker_contigs[marker_gene].append(contig_num)
                
                if marker_gene not in marker_contig_counts:
                    marker_contig_counts[marker_gene] = 1
                else:
                    marker_contig_counts[marker_gene] += 1

marker_frequencies = {}

for marker in marker_contig_counts:
    
    if marker_contig_counts[marker] not in marker_frequencies:
        marker_frequencies[marker_contig_counts[marker]] = 1
    else:
        marker_frequencies[marker_contig_counts[marker]] += 1



# Get marker gene counts to make bins
#-----------------------------------------------------

print("\nDetermining which marker genes to consider...")

my_gene_counts = []

vals = list(marker_contig_counts.values())

data = np.array(vals)

values, counts = np.unique(data, return_counts=True)

max_count_gene = values[np.where(counts == np.amax(counts))]

print("max_count_gene:", max_count_gene)

min_count_index = -1
max_count_index = -1

my_gene_counts.append(max_count_gene[0])

max_count_gene_index = np.where(values == max_count_gene[0])[0][0]

print(counts[max_count_gene_index]*seed_mg_threshold)

for i in range(max_count_gene_index+1, len(values)):
    if counts[max_count_gene_index]*seed_mg_threshold <= counts[i]:
        max_count_index = i

for i in range(max_count_gene_index, 0, -1):
    if counts[max_count_gene_index]*seed_mg_threshold <= counts[i-1]:
        min_count_index = i-1
        
for i in range(min_count_index, max_count_index+1):
    if values[i] not in my_gene_counts:
        my_gene_counts.append(values[i])

my_gene_counts.sort(reverse=True)


# Get contigs containing each marker gene for each iteration
seed_iter = {}

n=0
for i in range(len(my_gene_counts)):
    for item in marker_contig_counts:
        if marker_contig_counts[item] == my_gene_counts[i]:
            seed_iter[n] = marker_contigs[item]
            n+=1


# Initialise bins
#-----------------------------------------------------

bins = {}
bin_of_contig = {}
binned_contigs = []

binned_contigs_with_markers = []

print("\nInitialising bins...")

print("\nseed_iter[0]", seed_iter[0])

print()

for i in range(len(seed_iter[0])):
    
    binned_contigs_with_markers.append(seed_iter[0][i])

    contig_num = int(graph_to_contig_map_rev[seed_iter[0][i]])
    
    bins[i] = [contigs_map_rev[contig_num]]
    bin_of_contig[contigs_map_rev[contig_num]] = i
    binned_contigs.append(contigs_map_rev[contig_num])
  
print(bins)


# Assign to bins
#-----------------------------------------------------

def normpdf(x, mean, sd):
    var = float(sd)**2
    denom = sd*(2*math.pi)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom


edge_weights_per_iteration = {}

print("\nMatching and assigning contigs with marker genes to bins...")
        
for i in range(len(seed_iter)):
    
    print("Iteration", i, ":", len(seed_iter[i]), "contigs")
    
    if i>0:
        
        B = nx.Graph()
                
        common = set(binned_contigs_with_markers).intersection(set(seed_iter[i]))
        
        to_bin = list(set(seed_iter[i]) - common)
        
        my_seeds = []
        
        n_bins = len(bins)
        
        bottom_nodes = []

        for n in range(n_bins):
            contigid = bins[n][0]
            if contigid not in bottom_nodes:
                bottom_nodes.append(contigid)
               
        top_nodes = []
        edges = []
        
        if len(to_bin)!=0:

            n_bins = len(bins)
        
            for contig in to_bin:

                contig_num = int(graph_to_contig_map_rev[contig])

                contigid = contigs_map_rev[contig_num]

                if contigid not in top_nodes:
                    top_nodes.append(contigid)

                my_seeds.append(contigid)

                for b in range(n_bins):

                    log_prob_sum = 0
                    n_contigs = len(bins[b])

                    for j in range(n_contigs):

                        if contigid in tetramer_profiles and bins[b][j] in tetramer_profiles:

                            tetramer_dist = distance.euclidean(tetramer_profiles[contigid], tetramer_profiles[bins[b][j]])
                            prob_comp = normpdf(tetramer_dist, mu_intra, sigma_intra)/(normpdf(tetramer_dist, mu_intra, sigma_intra)+normpdf(tetramer_dist, mu_inter, sigma_inter))
                            
                            if contigid not in coverages or bins[b][j] not in coverages:
                                prob_cov = 0
                            else:
                                prob_cov = poisson.pmf(coverages[contigid], coverages[bins[b][j]])

                            prob_product = prob_comp * prob_cov

                            log_prob = 0

                            if prob_product!=0.0:
                                log_prob = - (math.log(prob_comp, 10) + math.log(prob_cov, 10))
                            else:
                                log_prob = max_weight

                            log_prob_sum += log_prob

                    if log_prob_sum != float("inf"):
                        edges.append((bins[b][0], contigid, log_prob_sum/n_contigs))
                    else:
                        edges.append((bins[b][0], contigid, max_weight))

            B.add_nodes_from(top_nodes, bipartite=0)
            B.add_nodes_from(bottom_nodes, bipartite=1)

            edge_weights = {}


            # Add edges only between nodes of opposite node sets
            for edge in edges:
                edge_weights[(edge[0], edge[1])] = edge[2]
                B.add_edge(edge[0], edge[1], weight = edge[2])

            edge_weights_per_iteration[i] = edge_weights

            top_nodes = {n for n, d in B.nodes(data=True) if d['bipartite']==0}
            bottom_nodes = set(B) - top_nodes

            my_matching = nx.algorithms.bipartite.matching.minimum_weight_full_matching(B, top_nodes, "weight")

            not_binned = {}

            matched_in_seed_iter = []

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
                                               
                        if edge_weights[(l, my_matching[l])]<=alpha_intra and math.floor(avg_path_len)<=dist_intra:
                            if my_matching[l] in contigs_map:
                                if contigs_map[my_matching[l]] in graph_to_contig_map:
                                    bins[b].append(my_matching[l])
                                    bin_of_contig[my_matching[l]] = b
                                    matched_in_seed_iter.append(my_matching[l])
                                    binned_contigs_with_markers.append(graph_to_contig_map[contigs_map[my_matching[l]]])
                                    binned_contigs.append(my_matching[l])

                        else:
                            not_binned[my_matching[l]] = (l, b)

            for nb in not_binned:
                    
                if edge_weights_per_iteration[i][(not_binned[nb][0], nb)]>=alpha_inter:
                    
                    path_len_sum = 0
                    
                    for contig_in_bin in bins[not_binned[nb][1]]:
                        shortest_paths = assembly_graph.get_shortest_paths(nb, to=contig_in_bin)
                        
                        if len(shortest_paths) != 0:
                            path_len_sum += len(shortest_paths[0])
                    
                    avg_path_len = path_len_sum/len(bins[not_binned[nb][1]])
                    
                    if math.floor(avg_path_len) >= dist_inter:
                        print("Creating new bin...")
                        bins[n_bins]=[nb]
                        bin_of_contig[nb] = n_bins
                        binned_contigs.append(nb)
                        n_bins+=1
                        binned_contigs_with_markers.append(graph_to_contig_map[contigs_map[nb]])

print()

for b in bins:
    print(b, ":", bins[b])


# Get binned and unbinned contigs
#-----------------------------------------------------

binned_contigs = list(bin_of_contig.keys())
    
unbinned_contigs = list(set([x for x in range(node_count)]) - set(binned_contigs))

print("\nNumber of binned contigs:", len(binned_contigs))
print("Number of unbinned contigs:", len(unbinned_contigs))


# Get isolated vertices and components without labels
#-----------------------------------------------------

isolated=[]

for i in range(node_count):
    
    neighbours = assembly_graph.neighbors(i, mode=ALL)
    
    if len(neighbours)==0:
        isolated.append(i)


non_isolated = []

for i in range(node_count):
    
    if i not in non_isolated and i in binned_contigs:

        component = []
        component.append(i)
        length = len(component)
        neighbours = assembly_graph.neighbors(i, mode=ALL)

        for neighbor in neighbours:
            if neighbor not in component:
                component.append(neighbor)

        component = list(set(component))

        while length!= len(component):

            length = len(component)

            for j in component:

                neighbours = assembly_graph.neighbors(j, mode=ALL)

                for neighbor in neighbours:
                    if neighbor not in component:
                        component.append(neighbor)

        labelled = False
        for j in component:
            if j in binned_contigs:
                labelled = True
                break

        if labelled:
            for j in component:
                if j not in non_isolated:
                    non_isolated.append(j)

print("\nNumber of non-isolated contigs:", len(non_isolated))

non_isolated_unbinned = list(set(non_isolated).intersection(set(unbinned_contigs)))

print("Number of non-isolated unbinned contigs:", len(non_isolated_unbinned))


# The BFS function to search labelled nodes
#-----------------------------------------------------

def runBFS(node, threhold=depth):
    queue = []
    visited = set()
    queue.append(node)
    depth = {}
    
    depth[node] = 0
    
    labelled_nodes = set()
    
    while (len(queue) > 0):
        active_node = queue.pop(0)
        visited.add(active_node)
        
        if active_node in binned_contigs and len(visited) > 1:
            
            # Get the bin of the current contig
            contig_bin = bin_of_contig[active_node]

            if node in tetramer_profiles and active_node in tetramer_profiles:
            
                tetramer_dist = distance.euclidean(tetramer_profiles[node], tetramer_profiles[active_node])
                prob_comp = normpdf(tetramer_dist, mu_intra, sigma_intra)/(normpdf(tetramer_dist, mu_intra, sigma_intra)+normpdf(tetramer_dist, mu_inter, sigma_inter))

            else:
                prob_comp = 0                               
            
            if node not in coverages or active_node not in coverages:
                prob_cov = 0
            else:
                prob_cov = poisson.pmf(coverages[node], coverages[active_node])
            
            if contig_lengths[node] >= min_length and contig_lengths[active_node] >= min_length:
                if prob_cov!=0.0 and prob_comp!=0.0:
                    labelled_nodes.add((node, active_node, contig_bin, depth[active_node], -(math.log(prob_cov, 10)), -math.log(prob_comp, 10)))
                elif prob_cov==0.0 and prob_comp!=0.0:
                    labelled_nodes.add((node, active_node, contig_bin, depth[active_node], max_weight, -math.log(prob_comp, 10)))
                elif prob_cov!=0.0 and prob_comp==0.0:
                    labelled_nodes.add((node, active_node, contig_bin, depth[active_node], -math.log(prob_cov, 10), max_weight))
                else:
                    labelled_nodes.add((node, active_node, contig_bin, depth[active_node], max_weight, max_weight))
            else:
                if prob_cov!=0.0:
                    labelled_nodes.add((node, active_node, contig_bin, depth[active_node], -math.log(prob_cov, 10), max_weight))
                else:
                    labelled_nodes.add((node, active_node, contig_bin, depth[active_node], max_weight, max_weight))
                
        else:
            for neighbour in assembly_graph.neighbors(active_node, mode=ALL):
                if neighbour not in visited:
                    depth[neighbour] = depth[active_node] + 1
                    if depth[neighbour] > threhold:
                        continue
                    queue.append(neighbour)
                    
    return labelled_nodes


# Propagate labels to unlabelled vertices
#-----------------------------------------------------

print("\nPropagating labels to unlabelled vertices...")

# Initialise progress bar
pbar = tqdm(total=len(non_isolated_unbinned))

class DataWrap:
    def __init__(self, data):
        self.data = data
        
    def __lt__(self, other):
        return (self.data[3], self.data[4], self.data[5])  < (other.data[3], other.data[4], other.data[5]) 
    
contigs_to_bin = set()

for contig in binned_contigs:
    if contig in non_isolated:
        closest_neighbours = filter(lambda x: x not in binned_contigs, assembly_graph.neighbors(contig, mode=ALL))
        contigs_to_bin.update(closest_neighbours)


sorted_node_list = []
sorted_node_list_ = [list(runBFS(x, threhold=depth)) for x in contigs_to_bin]
sorted_node_list_ = [item for sublist in sorted_node_list_ for item in sublist]

for data in sorted_node_list_:
    heapObj = DataWrap(data)
    heapq.heappush(sorted_node_list, heapObj)


while sorted_node_list:
    best_choice = heapq.heappop(sorted_node_list)    
    to_bin, binned, bin_, dist, cov_diff, comp_diff = best_choice.data
    
    
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
            candidates = list(runBFS(n, threhold=depth))
            for c in candidates:
                heapq.heappush(sorted_node_list, DataWrap(c))

# Close progress bar
pbar.close()


print("\nRemaining number of unbinned contigs:", len(unbinned_contigs))
print("Total number of binned contigs:", len(binned_contigs))


# Bin remaining unbinned contigs
#-----------------------------------

print("\nPropagating labels to remaining unlabelled vertices...")

def assign(contig):

    if contig_lengths[contig] >= min_length:

        min_log_prob = max_weight
        min_log_prob_bin = -1

        for b in range(len(bins)):
            
            log_prob_final = 0
            log_prob_sum = 0
            n_contigs = 0

            for j in range(len(bins[b])):

                if contig_lengths[bins[b][j]] >= min_length:

                    tetramer_dist = distance.euclidean(tetramer_profiles[contig], tetramer_profiles[bins[b][j]])
                    prob_comp = normpdf(tetramer_dist, mu_intra, sigma_intra)/(normpdf(tetramer_dist, mu_intra, sigma_intra)+normpdf(tetramer_dist, mu_inter, sigma_inter))
                    
                    # if contig not in coverages and bins[b][j] in coverages:
                    #     prob_cov = poisson.pmf(0, coverages[bins[b][j]])
                    # elif contig in coverages and bins[b][j] not in coverages:
                    #     prob_cov = poisson.pmf(coverages[contig], 0)
                    if contig not in coverages or bins[b][j] not in coverages:
                        continue
                    else:
                        prob_cov = poisson.pmf(coverages[contig], coverages[bins[b][j]])

                    prob_product = prob_comp * prob_cov

                    log_prob = 0

                    if prob_product!=0.0:
                        log_prob = - (math.log(prob_comp, 10) + math.log(prob_cov, 10))
                    else:
                        log_prob = max_weight

                    log_prob_sum += log_prob
                    n_contigs += 1

            if log_prob_sum != float("inf") and n_contigs!=0:
                log_prob_final = log_prob_sum/n_contigs
            else:
                log_prob_final = max_weight

            if min_log_prob > log_prob_final:
                min_log_prob = log_prob_final
                min_log_prob_bin = b

        if min_log_prob_bin !=-1:
            return contig, min_log_prob_bin

    return None


with Pool(nthreads) as p:
    assigned = list(tqdm(p.imap(assign, unbinned_contigs), total=len(unbinned_contigs)))

put_to_bins = list(filter(lambda x: x is not None, assigned))

if len(put_to_bins) ==  0:
    print("No further contigs were binned")
else:
    print(str(len(put_to_bins)), "contigs were binned")

# Assign contigs to bins
for contig, contig_bin in put_to_bins:
    bins[contig_bin].append(contig)
    bin_of_contig[contig] = contig_bin
    binned_contigs.append(contig)
    unbinned_contigs.remove(contig)

print("\nRemaining number of unbinned contigs:", len(unbinned_contigs))
print("Total number of binned contigs:", len(binned_contigs))


# Get elapsed time
#-----------------------------------

# Determine elapsed time
elapsed_time = time.time() - start_time

# Print elapsed time for the process
print("\nElapsed time: ", elapsed_time, " seconds")

# Sort contigs in bins
for i in range(n_bins):
    bins[i].sort()


# Write result to output file
#-----------------------------------

output_bins = []

for contig in bin_of_contig:
    line = []
    line.append(graph_to_contig_map[contigs_map[contig]])
    line.append(bin_of_contig[contig]+1)
    output_bins.append(line)

output_file = output_path + prefix + 'metacoag_output.csv'

with open(output_file, mode='w') as output_file:
    output_writer = csv.writer(output_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    for row in output_bins:
        output_writer.writerow(row)

print("\nFinal binning results can be found at", output_file.name)


# Exit program
#-----------------------------------

print("\nThank you for using MetaCoAG!\n")