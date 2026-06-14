#!/usr/bin/env python3

import concurrent.futures
import logging
import math
import operator
import sys

import networkx as nx
import numpy as np
from scipy.spatial import distance
from scipy.special import gammaln

__author__ = "Vijini Mallawaarachchi and Yu Lin"
__copyright__ = "Copyright 2020, MetaCoAG Project"
__license__ = "GPL-3.0"
__version__ = "1.2.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Stable Release"


# Constants set from MaxBin 2.0
MU_INTRA, SIGMA_INTRA = 0, 0.01037897 / 2
MU_INTER, SIGMA_INTER = 0.0676654, 0.03419337
VERY_SMALL_DOUBLE = 1e-10
MAX_WEIGHT = sys.float_info.max

# create logger
logger = logging.getLogger(f"MetaCoaAG {__version__}")

_WORKER_BIN_TETRA_MAT = None
_WORKER_BIN_COV_MAT = None
_WORKER_W_INTRA = None


def normpdf(x, mean, sd):
    # Vectorised via numpy; scalar inputs also work.
    x = np.asarray(x, dtype=float)
    return np.exp(-0.5 * ((x - mean) / sd) ** 2) / (sd * np.sqrt(2.0 * np.pi))


def get_tetramer_distance(seq1, seq2):
    return distance.euclidean(seq1, seq2)


def get_coverage_distance(cov1, cov2):
    return distance.euclidean(cov1, cov2)


def get_comp_probability(tetramer_dist):
    gaus_intra = normpdf(tetramer_dist, MU_INTRA, SIGMA_INTRA)
    gaus_inter = normpdf(tetramer_dist, MU_INTER, SIGMA_INTER)
    return float(gaus_intra / (gaus_intra + gaus_inter))


def get_cov_probability(cov1, cov2):
    # Vectorised Poisson PMF computation.
    # Adapted from http://www.masaers.com/2013/10/08/Implementing-Poisson-pmf.html
    c1 = np.asarray(cov1, dtype=float)
    c2 = np.asarray(cov2, dtype=float)
    # Guard against log(0): replace zeros with a tiny positive value
    safe_c1 = np.where(c1 > 0, c1, VERY_SMALL_DOUBLE)
    safe_c2 = np.where(c2 > 0, c2, VERY_SMALL_DOUBLE)
    log_pmf_1 = c1 * np.log(safe_c2) - gammaln(c1 + 1.0) - c2
    log_pmf_2 = c2 * np.log(safe_c1) - gammaln(c2 + 1.0) - c1
    pmf_1 = np.maximum(np.exp(log_pmf_1), VERY_SMALL_DOUBLE)
    pmf_2 = np.maximum(np.exp(log_pmf_2), VERY_SMALL_DOUBLE)
    return float(min(np.prod(pmf_1), np.prod(pmf_2)))


def _build_bin_matrices(bins, n_bins, normalized_tetramer_profiles, coverages):
    """Stack bin members into numpy arrays for vectorised per-iteration scoring."""
    bin_tetra_mat = {}
    bin_cov_mat = {}
    for b in range(n_bins):
        members = np.asarray(bins[b], dtype=np.intp)
        bin_tetra_mat[b] = normalized_tetramer_profiles[members]
        bin_cov_mat[b] = coverages[members]
    return bin_tetra_mat, bin_cov_mat


def _compute_edge_weight_exact(tetra_contig, cov_contig, bin_tetra_mat, bin_cov_mat):
    """Vectorised, exact equivalent of the original per-member scoring loop.

    Computes mean(-log10(p_comp_j) - log10(p_cov_j)) across all N bin members
    using numpy/scipy in C — identical numerical result to the Python for-j loop,
    including the float-overflow -> MAX_WEIGHT behaviour.
    """
    # All N tetramer distances in one cdist sweep
    dists = distance.cdist([tetra_contig], bin_tetra_mat, "euclidean")[0]  # (N,)

    # Composition probabilities for all N members simultaneously
    gi = np.exp(-0.5 * (dists / SIGMA_INTRA) ** 2) / (SIGMA_INTRA * np.sqrt(2.0 * np.pi))
    ge = np.exp(-0.5 * ((dists - MU_INTER) / SIGMA_INTER) ** 2) / (SIGMA_INTER * np.sqrt(2.0 * np.pi))
    prob_comp_vec = gi / (gi + ge)  # (N,)

    # Coverage probabilities: vectorised Poisson PMF over all N members
    c1 = np.asarray(cov_contig, dtype=float)              # (S,)
    mat_c = bin_cov_mat                                   # (N, S)
    s1 = np.where(c1 > 0, c1, VERY_SMALL_DOUBLE)
    s2 = np.where(mat_c > 0, mat_c, VERY_SMALL_DOUBLE)
    lp1 = c1 * np.log(s2) - gammaln(c1 + 1.0) - mat_c   # (N, S)
    lp2 = mat_c * np.log(s1) - gammaln(mat_c + 1.0) - c1 # (N, S)
    prod1 = np.prod(np.maximum(np.exp(lp1), VERY_SMALL_DOUBLE), axis=1)  # (N,)
    prod2 = np.prod(np.maximum(np.exp(lp2), VERY_SMALL_DOUBLE), axis=1)  # (N,)
    prob_cov_vec = np.minimum(prod1, prod2)               # (N,)

    # Per-member log probabilities — same formula as the original scalar path
    prob_product_vec = prob_comp_vec * prob_cov_vec
    mask = prob_product_vec > 0.0
    log_probs = np.where(
        mask,
        -(np.log10(np.where(mask, prob_comp_vec, 1.0)) +
          np.log10(np.where(mask, prob_cov_vec, 1.0))),
        MAX_WEIGHT,
    )  # (N,)

    # Reproduce original overflow check: any MAX_WEIGHT entry pushes the sum
    # to inf, causing the same MAX_WEIGHT result as the original loop.
    log_prob_sum = float(np.sum(log_probs))
    if math.isinf(log_prob_sum):
        return MAX_WEIGHT
    return log_prob_sum / len(bin_tetra_mat)


def match_contigs(
    smg_iteration,
    bins,
    n_bins,
    bin_of_contig,
    binned_contigs_with_markers,
    bin_markers,
    contig_markers,
    contig_lengths,
    contig_names,
    normalized_tetramer_profiles,
    coverages,
    assembly_graph,
    w_intra,
    w_inter,
    d_limit,
):
    edge_weights_per_iteration = {}

    smg_iterations = len(smg_iteration)

    for i in range(smg_iterations):
        logger.debug(
            f"Iteration {i}: {len(smg_iteration[i])} contig(s) with seed marker genes"
        )

        if i > 0:
            B = nx.Graph()

            common = set(binned_contigs_with_markers).intersection(
                set(smg_iteration[i])
            )
            to_bin = list(set(smg_iteration[i]) - common)
            logger.debug(f"{len(to_bin)} contig(s) to bin in the iteration")
            n_bins = len(bins)
            bottom_nodes = []

            for n in range(n_bins):
                contigid = bins[n][0]
                if contigid not in bottom_nodes:
                    bottom_nodes.append(contigid)

            top_nodes = []
            edges = []

            binned_count = 0

            # Build per-bin member matrices once per iteration so all members
            # assigned in previous iterations are included. Rebuilt next
            # iteration automatically.
            bin_tetra_mat, bin_cov_mat = _build_bin_matrices(
                bins, n_bins, normalized_tetramer_profiles, coverages
            )

            if len(to_bin) != 0:
                for contig in to_bin:
                    contigid = contig

                    if contigid not in top_nodes:
                        top_nodes.append(contigid)

                    tetra_contig = normalized_tetramer_profiles[contigid]
                    cov_contig = coverages[contigid]

                    for b in range(n_bins):
                        edge_weight = _compute_edge_weight_exact(
                            tetra_contig,
                            cov_contig,
                            bin_tetra_mat[b],
                            bin_cov_mat[b],
                        )
                        edges.append((bins[b][0], contigid, edge_weight))

                B.add_nodes_from(top_nodes, bipartite=0)
                B.add_nodes_from(bottom_nodes, bipartite=1)

                edge_weights = {}

                # Add edges only between nodes of opposite node sets
                for edge in edges:
                    edge_weights[(edge[0], edge[1])] = edge[2]
                    B.add_edge(edge[0], edge[1], weight=edge[2])

                edge_weights_per_iteration[i] = edge_weights

                top_nodes = {n for n, d in B.nodes(data=True) if d["bipartite"] == 0}
                bottom_nodes = set(B) - top_nodes

                if len(top_nodes) > 0:
                    my_matching = (
                        nx.algorithms.bipartite.matching.minimum_weight_full_matching(
                            B, top_nodes, "weight"
                        )
                    )

                    not_binned = {}

                    for l in my_matching:
                        if l in bin_of_contig:
                            b = bin_of_contig[l]

                            if (
                                my_matching[l] not in bins[b]
                                and (l, my_matching[l]) in edge_weights
                            ):
                                # Batch all targets in the bin into a single
                                # igraph distances() call (one BFS sweep).
                                all_paths = assembly_graph.distances(
                                    my_matching[l], target=bins[b]
                                )
                                # distances() returns a 2-D list; row 0 for our source
                                path_len_sum = sum(
                                    d for d in all_paths[0] if d != float("inf")
                                )

                                avg_path_len = math.floor(path_len_sum / len(bins[b]))

                                # logger.debug("To assign contig " + contig_names[my_matching[l]] + " to bin "+str(
                                #     b+1) + " based on contig " + str(l) + " weight="+str(edge_weights[(l, my_matching[l])]))

                                if (
                                    edge_weights[(l, my_matching[l])] <= w_intra
                                    and avg_path_len <= d_limit
                                ):
                                    can_assign = False

                                    common_mgs = (
                                        bin_markers[b]
                                        & contig_markers[my_matching[l]]
                                    )

                                    if len(common_mgs) == 0:
                                        can_assign = True
                                    # else:
                                    #     neighbours = assembly_graph.neighbors(my_matching[l])

                                    #     # if len(neighbours) == 0 and len(common_mgs) <= 1:
                                    #     if len(neighbours) == 0 and len(common_mgs) <= 1 and contig_lengths[my_matching[l]] > 2000:
                                    #         can_assign = True

                                    if can_assign:
                                        bins[b].append(my_matching[l])
                                        bin_of_contig[my_matching[l]] = b
                                        binned_contigs_with_markers.append(
                                            my_matching[l]
                                        )
                                        binned_count += 1

                                        # logger.debug("Assigning contig " + contig_names[my_matching[l]] + " to bin "+str(
                                        #     b+1) + " based on contig " + str(l) + " weight="+str(edge_weights[(l, my_matching[l])]))

                                        bin_markers[b].update(
                                            contig_markers[my_matching[l]]
                                        )

                                    else:
                                        not_binned[my_matching[l]] = (l, b)

                                else:
                                    not_binned[my_matching[l]] = (l, b)

                    longest_nb_contig = -1
                    longest_nb_contig_length = -1
                    longest_nb_contig_mg_count = -1

                    for nb in not_binned:
                        if (
                            edge_weights_per_iteration[i][(not_binned[nb][0], nb)]
                            > w_inter
                        ):
                            if longest_nb_contig_mg_count < len(
                                contig_markers[not_binned[nb][0]]
                            ):
                                longest_nb_contig = nb
                                longest_nb_contig_mg_count = len(
                                    contig_markers[not_binned[nb][0]]
                                )
                                longest_nb_contig_length = contig_lengths[
                                    not_binned[nb][0]
                                ]

                            elif longest_nb_contig_mg_count == len(
                                contig_markers[not_binned[nb][0]]
                            ):
                                if (
                                    longest_nb_contig_length
                                    < contig_lengths[not_binned[nb][0]]
                                ):
                                    longest_nb_contig = nb
                                    longest_nb_contig_mg_count = len(
                                        contig_markers[not_binned[nb][0]]
                                    )
                                    longest_nb_contig_length = contig_lengths[
                                        not_binned[nb][0]
                                    ]

                    if longest_nb_contig != -1:
                        target_bin = not_binned[longest_nb_contig][1]
                        all_paths = assembly_graph.distances(
                            longest_nb_contig, target=bins[target_bin]
                        )
                        path_len_sum = sum(
                            d for d in all_paths[0] if d != float("inf")
                        )

                        avg_path_len = path_len_sum / len(bins[target_bin])

                        if math.floor(avg_path_len) >= d_limit or path_len_sum == 0:
                            logger.debug("Creating new bin...")
                            logger.debug(
                                "New bin has contig "
                                + str(longest_nb_contig)
                                + " to bin "
                                + str(n_bins + 1)
                                + " weight="
                                + str(
                                    edge_weights_per_iteration[i][
                                        (
                                            not_binned[longest_nb_contig][0],
                                            longest_nb_contig,
                                        )
                                    ]
                                )
                            )
                            bins[n_bins] = [longest_nb_contig]
                            bin_of_contig[longest_nb_contig] = n_bins
                            binned_count += 1

                            bin_markers[n_bins] = contig_markers[
                                longest_nb_contig
                            ].copy()
                            n_bins += 1
                            binned_contigs_with_markers.append(longest_nb_contig)

            logger.debug(f"{binned_count} contig(s) binned in the iteration")

    if len(smg_iteration) > 0:
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

    return bins, bin_of_contig, n_bins, bin_markers, binned_contigs_with_markers


def _score_contig_against_bins_exact(
    contigid,
    possible_bins,
    tetra_contig,
    cov_contig,
    bin_tetra_mat,
    bin_cov_mat,
    w_intra,
):
    """Score one contig against its candidate bins."""
    bin_weights = [
        _compute_edge_weight_exact(
            tetra_contig,
            cov_contig,
            bin_tetra_mat[b],
            bin_cov_mat[b],
        )
        for b in possible_bins
    ]

    min_b_index, min_b_value = min(
        enumerate(bin_weights), key=operator.itemgetter(1)
    )

    if min_b_value <= w_intra:
        return contigid, possible_bins[min_b_index], min_b_value
    return None


def _init_contig_scoring_worker(bin_tetra_mat, bin_cov_mat, w_intra):
    """Install read-only scoring data once in each worker process."""
    global _WORKER_BIN_TETRA_MAT, _WORKER_BIN_COV_MAT, _WORKER_W_INTRA

    _WORKER_BIN_TETRA_MAT = bin_tetra_mat
    _WORKER_BIN_COV_MAT = bin_cov_mat
    _WORKER_W_INTRA = w_intra


def _score_contig_against_bins_worker(args):
    """Process-pool wrapper using matrices initialized once per worker."""
    contigid, possible_bins, tetra_contig, cov_contig = args
    return _score_contig_against_bins_exact(
        contigid=contigid,
        possible_bins=possible_bins,
        tetra_contig=tetra_contig,
        cov_contig=cov_contig,
        bin_tetra_mat=_WORKER_BIN_TETRA_MAT,
        bin_cov_mat=_WORKER_BIN_COV_MAT,
        w_intra=_WORKER_W_INTRA,
    )


def further_match_contigs(
    unbinned_mg_contigs,
    min_length,
    bins,
    n_bins,
    bin_of_contig,
    binned_contigs_with_markers,
    bin_markers,
    contig_markers,
    normalized_tetramer_profiles,
    coverages,
    w_intra,
    nthreads=1,
):
    # Build per-bin member matrices once. Parallel workers receive these once
    # during initialization instead of with every contig-scoring task.
    bin_tetra_mat, bin_cov_mat = _build_bin_matrices(
        bins, len(bins), normalized_tetramer_profiles, coverages
    )

    # Build work items for the parallel scoring phase.
    work_items = []
    for contig in unbinned_mg_contigs:
        if contig[1] < min_length:
            continue
        contigid = contig[0]
        contig_mg_set = contig_markers[contigid]
        possible_bins = [
            b for b in bin_markers if bin_markers[b].isdisjoint(contig_mg_set)
        ]
        if not possible_bins:
            continue
        work_items.append((
            contigid,
            possible_bins,
            normalized_tetramer_profiles[contigid],
            coverages[contigid],
        ))

    if nthreads == 1:
        results = [
            _score_contig_against_bins_exact(
                contigid=contigid,
                possible_bins=possible_bins,
                tetra_contig=tetra_contig,
                cov_contig=cov_contig,
                bin_tetra_mat=bin_tetra_mat,
                bin_cov_mat=bin_cov_mat,
                w_intra=w_intra,
            )
            for contigid, possible_bins, tetra_contig, cov_contig in work_items
        ]
    elif work_items:
        chunksize = max(1, math.ceil(len(work_items) / (nthreads * 4)))
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=nthreads,
            initializer=_init_contig_scoring_worker,
            initargs=(bin_tetra_mat, bin_cov_mat, w_intra),
        ) as executor:
            results = list(
                executor.map(
                    _score_contig_against_bins_worker,
                    work_items,
                    chunksize=chunksize,
                )
            )
    else:
        results = []

    for result in results:
        if result is None:
            continue
        contigid, best_bin, _ = result
        # Guard: contig may appear in multiple work items
        if contigid in bin_of_contig:
            continue
        bins[best_bin].append(contigid)
        bin_of_contig[contigid] = best_bin
        binned_contigs_with_markers.append(contigid)
        bin_markers[best_bin].update(contig_markers[contigid])

    return bins, bin_of_contig, n_bins, bin_markers, binned_contigs_with_markers
