#!/usr/bin/env python3

import sys
import math
import operator
import numpy as np

from metacoag_utils import matching_utils

MAX_WEIGHT = sys.float_info.max


class DataWrap:
    def __init__(self, data):
        self.data = data

    def __lt__(self, other):
        # return (self.data[3], self.data[4], self.data[5])  < (other.data[3], other.data[4], other.data[5])
        return (self.data[3], self.data[4]) < (other.data[3], other.data[4])


# The BFS function to search labelled nodes
def runBFSLong(
        node, threhold, min_length, binned_contigs, bin_of_contig,
        assembly_graph, tetramer_profiles, coverages, contig_lengths):

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

            tetramer_dist = matching_utils.get_tetramer_distance(
                tetramer_profiles[node], tetramer_profiles[active_node])

            prob_comp = matching_utils.get_comp_probability(tetramer_dist)
            prob_cov = matching_utils.get_cov_probability(
                coverages[node], coverages[active_node])

            prob_product = prob_comp * prob_cov

            log_prob = 0

            if prob_product > 0.0:
                log_prob = - (math.log(prob_comp, 10) + math.log(prob_cov, 10))
            else:
                log_prob = MAX_WEIGHT

            labelled_nodes.add((node, active_node, contig_bin, depth[active_node], log_prob))

        else:
            for neighbour in assembly_graph.neighbors(active_node, mode="ALL"):
                if neighbour not in visited:
                    depth[neighbour] = depth[active_node] + 1
                    if depth[neighbour] > threhold:
                        continue
                    queue.append(neighbour)

    return labelled_nodes


# The BFS function to search labelled nodes
# def runBFS(
#         node, threhold, binned_contigs, bin_of_contig,
#         assembly_graph, coverages):

#     queue = []
#     visited = set()
#     queue.append(node)
#     depth = {}

#     depth[node] = 0

#     labelled_nodes = set()

#     while (len(queue) > 0):
#         active_node = queue.pop(0)
#         visited.add(active_node)

#         if active_node in binned_contigs and len(visited) > 1:

#             # Get the bin of the current contig
#             contig_bin = bin_of_contig[active_node]

#             dist = np.linalg.norm(
#                 np. array(coverages[node]) - np. array(coverages[active_node]))

#             labelled_nodes.add(
#                 (node, active_node, contig_bin, depth[active_node], dist))

#         else:
#             for neighbour in assembly_graph.neighbors(active_node, mode="ALL"):
#                 if neighbour not in visited:
#                     depth[neighbour] = depth[active_node] + 1
#                     if depth[neighbour] > threhold:
#                         continue
#                     queue.append(neighbour)

#     return labelled_nodes


def getClosestLongVertices(graph, node, binned_contigs, contig_lengths, min_length):

    queu_l = [graph.neighbors(node, mode='ALL')]
    visited_l = [node]
    unlabelled = []

    while len(queu_l) > 0:
        active_level = queu_l.pop(0)
        is_finish = False
        visited_l += active_level

        for n in active_level:
            if contig_lengths[n] >= min_length and n not in binned_contigs:
                is_finish = True
                unlabelled.append(n)
        if is_finish:
            return unlabelled
        else:
            temp = []
            for n in active_level:
                temp += graph.neighbors(n, mode='ALL')
                temp = list(set(temp))
            temp2 = []

            for n in temp:
                if n not in visited_l:
                    temp2.append(n)
            if len(temp2) > 0:
                queu_l.append(temp2)
    return unlabelled


def assignLong(
        contigid, coverages, normalized_tetramer_profiles,
        bins, contig_lengths, seed_iters):

    bin_weights = []

    for b in bins:

        log_prob_sum = 0

        n_contigs = 0

        if len(bins[b]) > seed_iters:
            n_contigs = seed_iters
        else:
            n_contigs = len(bins[b])

        for j in range(len(bins[b])):

            tetramer_dist = matching_utils.get_tetramer_distance(normalized_tetramer_profiles[contigid],
                                                                normalized_tetramer_profiles[bins[b][j]])
            prob_comp = matching_utils.get_comp_probability(tetramer_dist)
            prob_cov = matching_utils.get_cov_probability(
                coverages[contigid], coverages[bins[b][j]])

            prob_product = prob_comp * prob_cov

            if prob_product > 0.0:
                log_prob_sum += - (math.log(prob_comp, 10) + math.log(prob_cov, 10))
                n_contigs += 1
            else:
                log_prob_sum = MAX_WEIGHT

        if log_prob_sum != float("inf") and n_contigs!=0:
            bin_weights.append(log_prob_sum/n_contigs)
        else:
            bin_weights.append(MAX_WEIGHT)

    min_b_index = -1

    min_b_index, min_b_value = min(
        enumerate(bin_weights), key=operator.itemgetter(1))

    if min_b_index != -1 and min_b_value != MAX_WEIGHT:
        return contigid, min_b_index, min_b_value

    return None
