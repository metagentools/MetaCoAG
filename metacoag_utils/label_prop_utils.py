#!/usr/bin/env python3

import heapq
import logging
import math
import sys

from metacoag_utils import matching_utils

MAX_WEIGHT = sys.float_info.max

# create logger
logger = logging.getLogger("MetaCoaAG 1.1.2")


class DataWrap:
    def __init__(self, data):
        self.data = data

    def __lt__(self, other):
        return (self.data[3], self.data[4]) < (other.data[3], other.data[4])


def run_bfs_long(
    node,
    threhold,
    binned_contigs,
    bin_of_contig,
    bins,
    smg_bin_counts,
    assembly_graph,
    normalized_tetramer_profiles,
    coverages,
):
    # Search labelled long contigs using BFS

    queue = []
    visited = set()
    queue.append(node)
    depth = {}

    depth[node] = 0

    labelled_nodes = set()

    while len(queue) > 0:
        active_node = queue.pop(0)
        visited.add(active_node)

        if active_node in binned_contigs and len(visited) > 1:
            # Get the bin of the current contig
            contig_bin = bin_of_contig[active_node]

            bin_log_prob = 0

            log_prob_sum = 0

            n_contigs = smg_bin_counts[contig_bin]
            bin_n_contigs = 0

            for j in range(n_contigs):
                tetramer_dist = matching_utils.get_tetramer_distance(
                    normalized_tetramer_profiles[node],
                    normalized_tetramer_profiles[bins[contig_bin][j]],
                )
                prob_comp = matching_utils.get_comp_probability(tetramer_dist)
                prob_cov = matching_utils.get_cov_probability(
                    coverages[node], coverages[bins[contig_bin][j]]
                )

                prob_product = prob_comp * prob_cov

                log_prob = 0

                if prob_product > 0.0:
                    log_prob = -(math.log(prob_comp, 10) + math.log(prob_cov, 10))
                    bin_n_contigs += 1
                else:
                    log_prob = MAX_WEIGHT

                log_prob_sum += log_prob

            if log_prob_sum != float("inf") and bin_n_contigs != 0:
                bin_log_prob = log_prob_sum / bin_n_contigs
            else:
                bin_log_prob = MAX_WEIGHT

            labelled_nodes.add(
                (node, active_node, contig_bin, depth[active_node], bin_log_prob)
            )

        else:
            for neighbour in assembly_graph.neighbors(active_node, mode="ALL"):
                if neighbour not in visited:
                    depth[neighbour] = depth[active_node] + 1
                    if depth[neighbour] > threhold:
                        continue
                    queue.append(neighbour)

    return labelled_nodes


def run_bfs_short(
    node, threhold, binned_contigs, bin_of_contig, assembly_graph, coverages
):
    # Search labelled contigs using BFS

    queue = []
    visited = set()
    queue.append(node)
    depth = {}

    depth[node] = 0

    labelled_nodes = set()

    while len(queue) > 0:
        active_node = queue.pop(0)
        visited.add(active_node)

        if active_node in binned_contigs and len(visited) > 1:
            # Get the bin of the current contig
            contig_bin = bin_of_contig[active_node]

            cov_dist = matching_utils.get_coverage_distance(
                coverages[active_node], coverages[node]
            )

            labelled_nodes.add(
                (node, active_node, contig_bin, depth[active_node], cov_dist)
            )

        else:
            for neighbour in assembly_graph.neighbors(active_node, mode="ALL"):
                if neighbour not in visited:
                    depth[neighbour] = depth[active_node] + 1
                    if depth[neighbour] > threhold:
                        continue
                    queue.append(neighbour)

    return labelled_nodes


def getClosestLongVertices(graph, node, binned_contigs, contig_lengths, min_length):
    queu_l = [graph.neighbors(node, mode="ALL")]
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
                temp += graph.neighbors(n, mode="ALL")
                temp = list(set(temp))
            temp2 = []

            for n in temp:
                if n not in visited_l:
                    temp2.append(n)
            if len(temp2) > 0:
                queu_l.append(temp2)
    return unlabelled


def label_prop(
    bin_of_contig,
    bins,
    contig_markers,
    bin_markers,
    binned_contigs_with_markers,
    smg_bin_counts,
    non_isolated,
    contig_lengths,
    min_length,
    assembly_graph,
    normalized_tetramer_profiles,
    coverages,
    depth,
    weight,
):
    contigs_to_bin = set()

    for contig in bin_of_contig:
        if contig in non_isolated and contig_lengths[contig] >= min_length:
            closest_neighbours = filter(
                lambda x: contig_lengths[x] >= min_length,
                getClosestLongVertices(
                    assembly_graph,
                    contig,
                    list(bin_of_contig.keys()),
                    contig_lengths,
                    min_length,
                ),
            )
            contigs_to_bin.update(closest_neighbours)

    sorted_node_list = []
    sorted_node_list_ = [
        list(
            run_bfs_long(
                x,
                depth,
                bin_of_contig.keys(),
                bin_of_contig,
                bins,
                smg_bin_counts,
                assembly_graph,
                normalized_tetramer_profiles,
                coverages,
            )
        )
        for x in contigs_to_bin
    ]
    sorted_node_list_ = [item for sublist in sorted_node_list_ for item in sublist]

    for data in sorted_node_list_:
        heapObj = DataWrap(data)
        heapq.heappush(sorted_node_list, heapObj)

    while sorted_node_list:
        best_choice = heapq.heappop(sorted_node_list)
        to_bin, binned, bin_, dist, cov_comp_diff = best_choice.data

        can_bin = False

        has_mg = False

        common_mgs = set()

        if to_bin in contig_markers:
            has_mg = True
            common_mgs = set(bin_markers[bin_]).intersection(
                set(contig_markers[to_bin])
            )

            if binned in contig_markers and dist == 1:
                neighbour_common_mgs = set(contig_markers[binned]).intersection(
                    set(contig_markers[to_bin])
                )

                if neighbour_common_mgs == common_mgs:
                    common_mgs = set()

        if to_bin not in bin_of_contig and cov_comp_diff < weight and dist <= depth:
            if len(common_mgs) == 0:
                can_bin = True
            elif len(common_mgs) <= 1 and contig_lengths[to_bin] > 100000:
                can_bin = True

        if can_bin:
            bins[bin_].append(to_bin)
            bin_of_contig[to_bin] = bin_

            if has_mg:
                binned_contigs_with_markers.append(to_bin)
                bin_markers[bin_] = list(
                    set(bin_markers[bin_] + contig_markers[to_bin])
                )

            # Discover to_bin's neighbours
            unbinned_neighbours = set(
                filter(
                    lambda x: contig_lengths[x] >= min_length,
                    getClosestLongVertices(
                        assembly_graph,
                        to_bin,
                        list(bin_of_contig.keys()),
                        contig_lengths,
                        min_length,
                    ),
                )
            )
            sorted_node_list = list(
                filter(lambda x: x.data[0] not in unbinned_neighbours, sorted_node_list)
            )
            heapq.heapify(sorted_node_list)

            for un in unbinned_neighbours:
                candidates = list(
                    run_bfs_long(
                        un,
                        depth,
                        list(bin_of_contig.keys()),
                        bin_of_contig,
                        bins,
                        smg_bin_counts,
                        assembly_graph,
                        normalized_tetramer_profiles,
                        coverages,
                    )
                )
                for c in candidates:
                    heapq.heappush(sorted_node_list, DataWrap(c))

    return bins, bin_of_contig, bin_markers, binned_contigs_with_markers


def assign_long(
    contigid,
    coverages,
    normalized_tetramer_profiles,
    bin_tetramer_profiles,
    bin_coverage_profiles,
):
    bin_weights = []

    # Get weight to each bin based on bin profiles
    for b in bin_tetramer_profiles:
        log_prob = 0

        tetramer_dist = matching_utils.get_tetramer_distance(
            normalized_tetramer_profiles[contigid], bin_tetramer_profiles[b]
        )
        prob_comp = matching_utils.get_comp_probability(tetramer_dist)
        prob_cov = matching_utils.get_cov_probability(
            coverages[contigid], bin_coverage_profiles[b]
        )

        prob_product = prob_comp * prob_cov

        if prob_product > 0.0:
            log_prob = -(math.log(prob_comp, 10) + math.log(prob_cov, 10))

        if log_prob != 0:
            bin_weights.append(log_prob)
        else:
            bin_weights.append(MAX_WEIGHT)

    # Get the bin with minimum weight
    min_index = [bin_weights.index(x) for x in sorted(bin_weights)[:1]][0]

    if bin_weights[min_index] != MAX_WEIGHT:
        return contigid, min_index, bin_weights[min_index]

    return None


def assign_to_bins(
    put_to_bins,
    bins,
    bin_of_contig,
    bin_markers,
    binned_contigs_with_markers,
    contig_markers,
    contig_lengths,
):
    for contig, min_index, bin_weight in put_to_bins:
        contig_bin = min_index

        if contig_bin is not None:
            can_bin = False

            has_mg = False

            common_mgs = []

            if contig in contig_markers:
                has_mg = True
                common_mgs = list(
                    set(bin_markers[contig_bin]).intersection(
                        set(contig_markers[contig])
                    )
                )

            if contig not in bin_of_contig and bin_weight != MAX_WEIGHT:
                if len(common_mgs) == 0:
                    can_bin = True
                elif len(common_mgs) <= 1 and contig_lengths[contig] > 100000:
                    can_bin = True

            if can_bin:
                bins[contig_bin].append(contig)
                bin_of_contig[contig] = contig_bin

                if has_mg:
                    binned_contigs_with_markers.append(contig)
                    bin_markers[contig_bin] = list(
                        set(bin_markers[contig_bin] + contig_markers[contig])
                    )

    return bins, bin_of_contig, bin_markers, binned_contigs_with_markers


def final_label_prop(
    bin_of_contig,
    bins,
    contig_markers,
    bin_markers,
    binned_contigs_with_markers,
    smg_bin_counts,
    contig_lengths,
    min_length,
    assembly_graph,
    normalized_tetramer_profiles,
    coverages,
    depth,
    weight,
):
    contigs_to_bin = set()

    for contig in bin_of_contig:
        if contig_lengths[contig] >= min_length:
            closest_neighbours = filter(
                lambda x: contig_lengths[x] >= min_length,
                getClosestLongVertices(
                    assembly_graph,
                    contig,
                    list(bin_of_contig.keys()),
                    contig_lengths,
                    min_length,
                ),
            )
            contigs_to_bin.update(closest_neighbours)

    sorted_node_list = []
    sorted_node_list_ = [
        list(
            run_bfs_long(
                x,
                depth,
                bin_of_contig.keys(),
                bin_of_contig,
                bins,
                smg_bin_counts,
                assembly_graph,
                normalized_tetramer_profiles,
                coverages,
            )
        )
        for x in contigs_to_bin
    ]
    sorted_node_list_ = [item for sublist in sorted_node_list_ for item in sublist]

    for data in sorted_node_list_:
        heapObj = DataWrap(data)
        heapq.heappush(sorted_node_list, heapObj)

    while sorted_node_list:
        best_choice = heapq.heappop(sorted_node_list)
        to_bin, binned, bin_, dist, cov_comp_diff = best_choice.data

        has_mg = False

        if to_bin in contig_markers:
            has_mg = True
            common_mgs = set(bin_markers[bin_]).intersection(
                set(contig_markers[to_bin])
            )

            if binned in contig_markers and dist == 1:
                neighbour_common_mgs = set(contig_markers[binned]).intersection(
                    set(contig_markers[to_bin])
                )

                if neighbour_common_mgs == common_mgs:
                    common_mgs = set()

        if to_bin not in bin_of_contig and cov_comp_diff != weight:
            bins[bin_].append(to_bin)
            bin_of_contig[to_bin] = bin_

            if has_mg:
                binned_contigs_with_markers.append(to_bin)
                bin_markers[bin_] = list(
                    set(bin_markers[bin_] + contig_markers[to_bin])
                )

            # Discover to_bin's neighbours
            unbinned_neighbours = set(
                filter(
                    lambda x: contig_lengths[x] >= min_length,
                    getClosestLongVertices(
                        assembly_graph,
                        to_bin,
                        list(bin_of_contig.keys()),
                        contig_lengths,
                        min_length,
                    ),
                )
            )
            sorted_node_list = list(
                filter(lambda x: x.data[0] not in unbinned_neighbours, sorted_node_list)
            )
            heapq.heapify(sorted_node_list)

            for un in unbinned_neighbours:
                candidates = list(
                    run_bfs_long(
                        un,
                        depth,
                        list(bin_of_contig.keys()),
                        bin_of_contig,
                        bins,
                        smg_bin_counts,
                        assembly_graph,
                        normalized_tetramer_profiles,
                        coverages,
                    )
                )
                for c in candidates:
                    heapq.heappush(sorted_node_list, DataWrap(c))

    return bins, bin_of_contig, bin_markers, binned_contigs_with_markers
