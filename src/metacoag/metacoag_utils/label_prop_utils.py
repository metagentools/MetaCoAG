#!/usr/bin/env python3

import concurrent.futures
import heapq
import logging
import math
import sys
from collections import deque

import numpy as np

from metacoag.metacoag_utils import matching_utils

MAX_WEIGHT = sys.float_info.max

__author__ = "Vijini Mallawaarachchi and Yu Lin"
__copyright__ = "Copyright 2020, MetaCoAG Project"
__license__ = "GPL-3.0"
__version__ = "1.2.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Stable Release"


# create logger
logger = logging.getLogger(f"MetaCoaAG {__version__}")


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
    bin_tetra_mat=None,
    bin_cov_mat=None,
):
    # Search labelled long contigs using BFS

    queue = deque()
    visited = set()
    queue.append(node)
    depth = {}

    depth[node] = 0

    labelled_nodes = set()

    while queue:
        active_node = queue.popleft()
        visited.add(active_node)

        if active_node in binned_contigs and len(visited) > 1:
            # Get the bin of the current contig
            contig_bin = bin_of_contig[active_node]

            if bin_tetra_mat is not None and contig_bin in bin_tetra_mat:
                # Vectorised path: one cdist call for all N seed members at once.
                # Mathematically identical to the scalar loop below —
                # same formula, same overflow-to-MAX_WEIGHT behaviour.
                bin_log_prob = matching_utils._compute_edge_weight_exact(
                    normalized_tetramer_profiles[node],
                    coverages[node],
                    bin_tetra_mat[contig_bin],
                    bin_cov_mat[contig_bin],
                )
            else:
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

    queue = deque()
    visited = set()
    queue.append(node)
    depth = {}

    depth[node] = 0

    labelled_nodes = set()

    while queue:
        active_node = queue.popleft()
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
    # binned_contigs must support O(1) membership tests (set or dict)
    queu_l = deque([graph.neighbors(node, mode="ALL")])
    visited_l = {node}
    unlabelled = []

    while queu_l:
        active_level = queu_l.popleft()
        is_finish = False
        visited_l.update(active_level)

        for n in active_level:
            if contig_lengths[n] >= min_length and n not in binned_contigs:
                is_finish = True
                unlabelled.append(n)
        if is_finish:
            return unlabelled
        else:
            temp = set()
            for n in active_level:
                temp.update(graph.neighbors(n, mode="ALL"))
            temp2 = [n for n in temp if n not in visited_l]
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
    nthreads=1,
):
    contigs_to_bin = set()

    # Use bin_of_contig directly (dict) for O(1) membership in getClosestLongVertices
    for contig in bin_of_contig:
        if contig in non_isolated and contig_lengths[contig] >= min_length:
            closest_neighbours = filter(
                lambda x: contig_lengths[x] >= min_length,
                getClosestLongVertices(
                    assembly_graph,
                    contig,
                    bin_of_contig,
                    contig_lengths,
                    min_length,
                ),
            )
            contigs_to_bin.update(closest_neighbours)

    sorted_node_list = []
    # Build seed-member matrices once. smg_bin_counts is frozen before label
    # propagation starts (computed from the initial seed bins) and never updated
    # as new contigs are added — so bins[b][:smg_bin_counts[b]] is stable
    # throughout the entire function, including the per-neighbour BFS calls
    # inside the assignment loop.
    _seed_tetra_mat = {}
    _seed_cov_mat = {}
    for _b in range(len(smg_bin_counts)):
        _n = smg_bin_counts[_b]
        _members = np.asarray(bins[_b][:_n], dtype=np.intp)
        _seed_tetra_mat[_b] = normalized_tetramer_profiles[_members]
        _seed_cov_mat[_b] = coverages[_members]

    # All BFS calls are independent (read-only data); run them in parallel.
    _binned_view = bin_of_contig.keys()
    def _bfs_long_worker_lp(x):
        return list(run_bfs_long(
            x, depth, _binned_view, bin_of_contig, bins, smg_bin_counts,
            assembly_graph, normalized_tetramer_profiles, coverages,
            bin_tetra_mat=_seed_tetra_mat, bin_cov_mat=_seed_cov_mat,
        ))
    with concurrent.futures.ThreadPoolExecutor(max_workers=nthreads) as pool:
        sorted_node_list_ = list(pool.map(_bfs_long_worker_lp, contigs_to_bin))
    sorted_node_list_ = [item for sublist in sorted_node_list_ for item in sublist]

    for data in sorted_node_list_:
        heapq.heappush(sorted_node_list, DataWrap(data))

    # Lazy-deletion set: contigs that already have a fresh BFS entry queued
    stale = set()

    while sorted_node_list:
        best_choice = heapq.heappop(sorted_node_list)
        to_bin, binned, bin_, dist, cov_comp_diff = best_choice.data

        # Skip stale entries whose neighbourhood has already been re-queued
        if to_bin in stale:
            continue

        can_bin = False

        has_mg = False

        common_mgs = set()

        if to_bin in contig_markers:
            has_mg = True
            common_mgs = bin_markers[bin_] & contig_markers[to_bin]

            if binned in contig_markers and dist == 1:
                neighbour_common_mgs = (
                    contig_markers[binned] & contig_markers[to_bin]
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
                bin_markers[bin_].update(contig_markers[to_bin])

            # Discover to_bin's neighbours; mark old entries stale instead of
            # rebuilding the heap, then push fresh BFS results.
            unbinned_neighbours = set(
                filter(
                    lambda x: contig_lengths[x] >= min_length,
                    getClosestLongVertices(
                        assembly_graph,
                        to_bin,
                        bin_of_contig,
                        contig_lengths,
                        min_length,
                    ),
                )
            )
            stale.update(unbinned_neighbours)

            for un in unbinned_neighbours:
                candidates = list(
                    run_bfs_long(
                        un,
                        depth,
                        bin_of_contig.keys(),
                        bin_of_contig,
                        bins,
                        smg_bin_counts,
                        assembly_graph,
                        normalized_tetramer_profiles,
                        coverages,
                        bin_tetra_mat=_seed_tetra_mat,
                        bin_cov_mat=_seed_cov_mat,
                    )
                )
                for c in candidates:
                    heapq.heappush(sorted_node_list, DataWrap(c))
                # Fresh entry is now queued; remove from stale so it can be processed
                stale.discard(un)

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
    min_index, min_weight = min(enumerate(bin_weights), key=lambda item: item[1])

    if min_weight != MAX_WEIGHT:
        return contigid, min_index, min_weight

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

            common_mgs = set()

            if contig in contig_markers:
                has_mg = True
                common_mgs = bin_markers[contig_bin] & contig_markers[contig]

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
                    bin_markers[contig_bin].update(contig_markers[contig])

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
    nthreads=1,
):
    contigs_to_bin = set()

    # Use bin_of_contig directly (dict) for O(1) membership in getClosestLongVertices
    for contig in bin_of_contig:
        if contig_lengths[contig] >= min_length:
            closest_neighbours = filter(
                lambda x: contig_lengths[x] >= min_length,
                getClosestLongVertices(
                    assembly_graph,
                    contig,
                    bin_of_contig,
                    contig_lengths,
                    min_length,
                ),
            )
            contigs_to_bin.update(closest_neighbours)

    sorted_node_list = []
    # Build seed-member matrices once. Same rationale as in label_prop: smg_bin_counts
    # is frozen so bins[b][:smg_bin_counts[b]] is stable throughout this function.
    _seed_tetra_mat_flp = {}
    _seed_cov_mat_flp = {}
    for _b in range(len(smg_bin_counts)):
        _n = smg_bin_counts[_b]
        _members = np.asarray(bins[_b][:_n], dtype=np.intp)
        _seed_tetra_mat_flp[_b] = normalized_tetramer_profiles[_members]
        _seed_cov_mat_flp[_b] = coverages[_members]

    # All BFS calls are independent (read-only data); run them in parallel.
    _binned_view_flp = bin_of_contig.keys()
    def _bfs_long_worker_flp(x):
        return list(run_bfs_long(
            x, depth, _binned_view_flp, bin_of_contig, bins, smg_bin_counts,
            assembly_graph, normalized_tetramer_profiles, coverages,
            bin_tetra_mat=_seed_tetra_mat_flp, bin_cov_mat=_seed_cov_mat_flp,
        ))
    with concurrent.futures.ThreadPoolExecutor(max_workers=nthreads) as pool:
        sorted_node_list_ = list(pool.map(_bfs_long_worker_flp, contigs_to_bin))
    sorted_node_list_ = [item for sublist in sorted_node_list_ for item in sublist]

    for data in sorted_node_list_:
        heapq.heappush(sorted_node_list, DataWrap(data))

    # Lazy-deletion set: contigs that already have a fresh BFS entry queued
    stale = set()

    while sorted_node_list:
        best_choice = heapq.heappop(sorted_node_list)
        to_bin, binned, bin_, dist, cov_comp_diff = best_choice.data

        # Skip stale entries whose neighbourhood has already been re-queued
        if to_bin in stale:
            continue

        has_mg = False

        if to_bin in contig_markers:
            has_mg = True
            common_mgs = bin_markers[bin_] & contig_markers[to_bin]

            if binned in contig_markers and dist == 1:
                neighbour_common_mgs = (
                    contig_markers[binned] & contig_markers[to_bin]
                )

                if neighbour_common_mgs == common_mgs:
                    common_mgs = set()

        if to_bin not in bin_of_contig and cov_comp_diff != weight:
            bins[bin_].append(to_bin)
            bin_of_contig[to_bin] = bin_

            if has_mg:
                binned_contigs_with_markers.append(to_bin)
                bin_markers[bin_].update(contig_markers[to_bin])

            # Discover to_bin's neighbours; mark old entries stale instead of
            # rebuilding the heap, then push fresh BFS results.
            unbinned_neighbours = set(
                filter(
                    lambda x: contig_lengths[x] >= min_length,
                    getClosestLongVertices(
                        assembly_graph,
                        to_bin,
                        bin_of_contig,
                        contig_lengths,
                        min_length,
                    ),
                )
            )
            stale.update(unbinned_neighbours)

            for un in unbinned_neighbours:
                candidates = list(
                    run_bfs_long(
                        un,
                        depth,
                        bin_of_contig.keys(),
                        bin_of_contig,
                        bins,
                        smg_bin_counts,
                        assembly_graph,
                        normalized_tetramer_profiles,
                        coverages,
                        bin_tetra_mat=_seed_tetra_mat_flp,
                        bin_cov_mat=_seed_cov_mat_flp,
                    )
                )
                for c in candidates:
                    heapq.heappush(sorted_node_list, DataWrap(c))
                # Fresh entry is now queued; remove from stale so it can be processed
                stale.discard(un)

    return bins, bin_of_contig, bin_markers, binned_contigs_with_markers
