#!/usr/bin/env python3

import sys
import math

from metacoag_utils import matching_utils

max_weight = sys.float_info.max


class DataWrap:
    def __init__(self, data):
        self.data = data
        
    def __lt__(self, other):
        # return (self.data[3], self.data[4], self.data[5])  < (other.data[3], other.data[4], other.data[5]) 
        return (self.data[3], self.data[4])  < (other.data[3], other.data[4]) 


# The BFS function to search labelled nodes
def runBFS(node, threhold, min_length, binned_contigs, bin_of_contig, assembly_graph, tetramer_profiles, coverages, contig_lengths):
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
            
            tetramer_dist = matching_utils.get_tetramer_distance(tetramer_profiles[node], tetramer_profiles[active_node])
                                
            prob_comp = matching_utils.get_comp_probability(tetramer_dist)
            prob_cov = matching_utils.get_cov_probability(coverages[node], coverages[active_node])

            if contig_lengths[node] >= min_length and contig_lengths[active_node] >= min_length:
                if prob_cov!=0.0 and prob_comp!=0.0:
                    labelled_nodes.add((node, active_node, contig_bin, depth[active_node], -(math.log(prob_cov, 10)+math.log(prob_comp, 10))))
                elif prob_cov==0.0 and prob_comp!=0.0:
                    labelled_nodes.add((node, active_node, contig_bin, depth[active_node], -math.log(prob_comp, 10)))
                elif prob_cov!=0.0 and prob_comp==0.0:
                    labelled_nodes.add((node, active_node, contig_bin, depth[active_node], -math.log(prob_cov, 10)))
                else:
                    labelled_nodes.add((node, active_node, contig_bin, depth[active_node], max_weight))
            else:
                if prob_cov!=0.0:
                    labelled_nodes.add((node, active_node, contig_bin, depth[active_node], -math.log(prob_cov, 10)))
                else:
                    labelled_nodes.add((node, active_node, contig_bin, depth[active_node], max_weight))
                
        else:
            for neighbour in assembly_graph.neighbors(active_node, mode="ALL"):
                if neighbour not in visited:
                    depth[neighbour] = depth[active_node] + 1
                    if depth[neighbour] > threhold:
                        continue
                    queue.append(neighbour)
                    
    return labelled_nodes



def assign(contig, min_length, contig_lengths, normalized_tetramer_profiles, bin_tetramer_profiles, n_bins):

    if contig_lengths[contig] >= min_length:

        min_tetramer_dist = max_weight
        min_tetramer_dist_bin = -1

        for b in range(n_bins):

            tetramer_dist = matching_utils.get_tetramer_distance(normalized_tetramer_profiles[contig], bin_tetramer_profiles[b])

            if min_tetramer_dist > tetramer_dist:
                min_tetramer_dist = tetramer_dist
                min_tetramer_dist_bin = b

        if min_tetramer_dist_bin !=-1:
            return contig, min_tetramer_dist_bin

    return None