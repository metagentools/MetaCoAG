#!/usr/bin/env python3

import sys
import math

from metacoag_utils import matching_utils

max_weight = sys.float_info.max


class DataWrap:
    def __init__(self, data):
        self.data = data
        
    def __lt__(self, other):
        return (self.data[3], self.data[4], self.data[5])  < (other.data[3], other.data[4], other.data[5]) 


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
            for neighbour in assembly_graph.neighbors(active_node, mode="ALL"):
                if neighbour not in visited:
                    depth[neighbour] = depth[active_node] + 1
                    if depth[neighbour] > threhold:
                        continue
                    queue.append(neighbour)
                    
    return labelled_nodes



def assign(contig, min_length, tetramer_profiles, coverages, contig_lengths, bins):

    if contig_lengths[contig] >= min_length:

        min_log_prob = max_weight
        min_log_prob_bin = -1

        for b in range(len(bins)):
            
            log_prob_final = 0
            log_prob_sum = 0
            n_contigs = 0

            for j in range(len(bins[b])):

                if contig_lengths[bins[b][j]] >= min_length:

                    tetramer_dist = matching_utils.get_tetramer_distance(tetramer_profiles[contig], tetramer_profiles[bins[b][j]])
                    prob_comp = matching_utils.get_comp_probability(tetramer_dist)
                    prob_cov = matching_utils.get_cov_probability(coverages[contig], coverages[bins[b][j]])

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