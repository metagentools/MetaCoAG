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

logger = logging.getLogger(f"MetaCoAG {__version__}")

_WORKER_BIN_TETRAMER_MATRICES = None
_WORKER_BIN_COVERAGE_MATRICES = None
_WORKER_W_INTRA = None


def normal_pdf(values, mean, standard_deviation):
    """Evaluate a normal probability density for scalar or array inputs."""
    values = np.asarray(values, dtype=float)
    return np.exp(-0.5 * ((values - mean) / standard_deviation) ** 2) / (
        standard_deviation * np.sqrt(2.0 * np.pi)
    )


normpdf = normal_pdf


def get_tetramer_distance(profile_a, profile_b):
    return distance.euclidean(profile_a, profile_b)


def get_coverage_distance(coverage_a, coverage_b):
    return distance.euclidean(coverage_a, coverage_b)


def get_comp_probability(tetramer_dist):
    intra_density = normal_pdf(tetramer_dist, MU_INTRA, SIGMA_INTRA)
    inter_density = normal_pdf(tetramer_dist, MU_INTER, SIGMA_INTER)
    return float(intra_density / (intra_density + inter_density))


def get_cov_probability(coverage_a, coverage_b):
    # Vectorised Poisson PMF computation.
    # Adapted from http://www.masaers.com/2013/10/08/Implementing-Poisson-pmf.html
    coverage_a = np.asarray(coverage_a, dtype=float)
    coverage_b = np.asarray(coverage_b, dtype=float)
    # Guard against log(0): replace zeros with a tiny positive value
    safe_coverage_a = np.where(coverage_a > 0, coverage_a, VERY_SMALL_DOUBLE)
    safe_coverage_b = np.where(coverage_b > 0, coverage_b, VERY_SMALL_DOUBLE)
    log_pmf_1 = (
        coverage_a * np.log(safe_coverage_b) - gammaln(coverage_a + 1.0) - coverage_b
    )
    log_pmf_2 = (
        coverage_b * np.log(safe_coverage_a) - gammaln(coverage_b + 1.0) - coverage_a
    )
    pmf_1 = np.maximum(np.exp(log_pmf_1), VERY_SMALL_DOUBLE)
    pmf_2 = np.maximum(np.exp(log_pmf_2), VERY_SMALL_DOUBLE)
    return float(min(np.prod(pmf_1), np.prod(pmf_2)))


def _build_bin_matrices(bins, n_bins, normalized_tetramer_profiles, coverages):
    """Stack bin members into numpy arrays for vectorised per-iteration scoring."""
    bin_tetramer_matrices = {}
    bin_coverage_matrices = {}
    for bin_id in range(n_bins):
        members = np.asarray(bins[bin_id], dtype=np.intp)
        bin_tetramer_matrices[bin_id] = normalized_tetramer_profiles[members]
        bin_coverage_matrices[bin_id] = coverages[members]
    return bin_tetramer_matrices, bin_coverage_matrices


def _compute_edge_weight_exact(
    contig_tetramers,
    contig_coverage,
    bin_tetramer_matrix,
    bin_coverage_matrix,
):
    """Vectorised, exact equivalent of the original per-member scoring loop.

    Computes mean(-log10(p_comp_j) - log10(p_cov_j)) across all N bin members
    using numpy/scipy in C — identical numerical result to the Python for-j loop,
    including the float-overflow -> MAX_WEIGHT behaviour.
    """
    # All N tetramer distances in one cdist sweep
    distances = distance.cdist([contig_tetramers], bin_tetramer_matrix, "euclidean")[0]

    # Composition probabilities for all N members simultaneously
    intra_densities = normal_pdf(distances, MU_INTRA, SIGMA_INTRA)
    inter_densities = normal_pdf(distances, MU_INTER, SIGMA_INTER)
    composition_probabilities = intra_densities / (intra_densities + inter_densities)

    # Coverage probabilities: vectorised Poisson PMF over all N members
    contig_coverage = np.asarray(contig_coverage, dtype=float)
    safe_contig_coverage = np.where(
        contig_coverage > 0, contig_coverage, VERY_SMALL_DOUBLE
    )
    safe_bin_coverages = np.where(
        bin_coverage_matrix > 0, bin_coverage_matrix, VERY_SMALL_DOUBLE
    )
    forward_log_pmf = (
        contig_coverage * np.log(safe_bin_coverages)
        - gammaln(contig_coverage + 1.0)
        - bin_coverage_matrix
    )
    reverse_log_pmf = (
        bin_coverage_matrix * np.log(safe_contig_coverage)
        - gammaln(bin_coverage_matrix + 1.0)
        - contig_coverage
    )
    forward_products = np.prod(
        np.maximum(np.exp(forward_log_pmf), VERY_SMALL_DOUBLE), axis=1
    )
    reverse_products = np.prod(
        np.maximum(np.exp(reverse_log_pmf), VERY_SMALL_DOUBLE), axis=1
    )
    coverage_probabilities = np.minimum(forward_products, reverse_products)

    # Per-member log probabilities — same formula as the original scalar path
    probability_products = composition_probabilities * coverage_probabilities
    valid_probabilities = probability_products > 0.0
    log_probs = np.where(
        valid_probabilities,
        -(
            np.log10(np.where(valid_probabilities, composition_probabilities, 1.0))
            + np.log10(np.where(valid_probabilities, coverage_probabilities, 1.0))
        ),
        MAX_WEIGHT,
    )

    # Reproduce original overflow check: any MAX_WEIGHT entry pushes the sum
    # to inf, causing the same MAX_WEIGHT result as the original loop.
    log_prob_sum = float(np.sum(log_probs))
    if math.isinf(log_prob_sum):
        return MAX_WEIGHT
    return log_prob_sum / len(bin_tetramer_matrix)


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
            bipartite_graph = nx.Graph()

            common = set(binned_contigs_with_markers).intersection(
                set(smg_iteration[i])
            )
            to_bin = list(set(smg_iteration[i]) - common)
            logger.debug(f"{len(to_bin)} contig(s) to bin in the iteration")
            n_bins = len(bins)
            bottom_nodes = []

            for bin_id in range(n_bins):
                seed_contig_id = bins[bin_id][0]
                if seed_contig_id not in bottom_nodes:
                    bottom_nodes.append(seed_contig_id)

            top_nodes = []
            edges = []

            binned_count = 0

            # Build per-bin member matrices once per iteration so all members
            # assigned in previous iterations are included. Rebuilt next
            # iteration automatically.
            bin_tetramer_matrices, bin_coverage_matrices = _build_bin_matrices(
                bins, n_bins, normalized_tetramer_profiles, coverages
            )

            if to_bin:
                for contig_id in to_bin:
                    if contig_id not in top_nodes:
                        top_nodes.append(contig_id)

                    contig_tetramers = normalized_tetramer_profiles[contig_id]
                    contig_coverage = coverages[contig_id]

                    for bin_id in range(n_bins):
                        edge_weight = _compute_edge_weight_exact(
                            contig_tetramers,
                            contig_coverage,
                            bin_tetramer_matrices[bin_id],
                            bin_coverage_matrices[bin_id],
                        )
                        edges.append((bins[bin_id][0], contig_id, edge_weight))

                bipartite_graph.add_nodes_from(top_nodes, bipartite=0)
                bipartite_graph.add_nodes_from(bottom_nodes, bipartite=1)

                edge_weights = {}

                # Add edges only between nodes of opposite node sets
                for edge in edges:
                    edge_weights[(edge[0], edge[1])] = edge[2]
                    bipartite_graph.add_edge(edge[0], edge[1], weight=edge[2])

                edge_weights_per_iteration[i] = edge_weights

                top_nodes = {
                    node
                    for node, data in bipartite_graph.nodes(data=True)
                    if data["bipartite"] == 0
                }
                bottom_nodes = set(bipartite_graph) - top_nodes

                if top_nodes:
                    matching = (
                        nx.algorithms.bipartite.matching.minimum_weight_full_matching(
                            bipartite_graph, top_nodes, "weight"
                        )
                    )

                    not_binned = {}

                    for source_contig in matching:
                        if source_contig in bin_of_contig:
                            bin_id = bin_of_contig[source_contig]
                            target_contig = matching[source_contig]

                            if (
                                target_contig not in bins[bin_id]
                                and (source_contig, target_contig) in edge_weights
                            ):
                                # Batch all targets in the bin into a single
                                # igraph distances() call (one BFS sweep).
                                all_paths = assembly_graph.distances(
                                    target_contig, target=bins[bin_id]
                                )
                                # distances() returns a 2-D list; row 0 for our source
                                path_len_sum = sum(
                                    d for d in all_paths[0] if d != float("inf")
                                )

                                avg_path_len = math.floor(
                                    path_len_sum / len(bins[bin_id])
                                )

                                if (
                                    edge_weights[(source_contig, target_contig)]
                                    <= w_intra
                                    and avg_path_len <= d_limit
                                ):
                                    can_assign = False

                                    common_mgs = (
                                        bin_markers[bin_id]
                                        & contig_markers[target_contig]
                                    )

                                    if len(common_mgs) == 0:
                                        can_assign = True

                                    if can_assign:
                                        bins[bin_id].append(target_contig)
                                        bin_of_contig[target_contig] = bin_id
                                        binned_contigs_with_markers.append(
                                            target_contig
                                        )
                                        binned_count += 1

                                        bin_markers[bin_id].update(
                                            contig_markers[target_contig]
                                        )

                                    else:
                                        not_binned[target_contig] = (
                                            source_contig,
                                            bin_id,
                                        )

                                else:
                                    not_binned[target_contig] = (
                                        source_contig,
                                        bin_id,
                                    )

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
                        path_len_sum = sum(d for d in all_paths[0] if d != float("inf"))

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

    return bins, bin_of_contig, n_bins, bin_markers, binned_contigs_with_markers


def _score_contig_against_bins_exact(
    contig_id,
    possible_bins,
    contig_tetramers,
    contig_coverage,
    bin_tetramer_matrices,
    bin_coverage_matrices,
    w_intra,
):
    """Score one contig against its candidate bins."""
    bin_weights = [
        _compute_edge_weight_exact(
            contig_tetramers,
            contig_coverage,
            bin_tetramer_matrices[bin_id],
            bin_coverage_matrices[bin_id],
        )
        for bin_id in possible_bins
    ]

    min_b_index, min_b_value = min(enumerate(bin_weights), key=operator.itemgetter(1))

    if min_b_value <= w_intra:
        return contig_id, possible_bins[min_b_index], min_b_value
    return None


def _init_contig_scoring_worker(bin_tetramer_matrices, bin_coverage_matrices, w_intra):
    """Install read-only scoring data once in each worker process."""
    global _WORKER_BIN_TETRAMER_MATRICES
    global _WORKER_BIN_COVERAGE_MATRICES
    global _WORKER_W_INTRA

    _WORKER_BIN_TETRAMER_MATRICES = bin_tetramer_matrices
    _WORKER_BIN_COVERAGE_MATRICES = bin_coverage_matrices
    _WORKER_W_INTRA = w_intra


def _score_contig_against_bins_worker(args):
    """Process-pool wrapper using matrices initialized once per worker."""
    contig_id, possible_bins, contig_tetramers, contig_coverage = args
    return _score_contig_against_bins_exact(
        contig_id=contig_id,
        possible_bins=possible_bins,
        contig_tetramers=contig_tetramers,
        contig_coverage=contig_coverage,
        bin_tetramer_matrices=_WORKER_BIN_TETRAMER_MATRICES,
        bin_coverage_matrices=_WORKER_BIN_COVERAGE_MATRICES,
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
    bin_tetramer_matrices, bin_coverage_matrices = _build_bin_matrices(
        bins, len(bins), normalized_tetramer_profiles, coverages
    )

    # Build work items for the parallel scoring phase.
    work_items = []
    for contig in unbinned_mg_contigs:
        if contig[1] < min_length:
            continue
        contig_id = contig[0]
        contig_marker_set = contig_markers[contig_id]
        possible_bins = [
            bin_id
            for bin_id in bin_markers
            if bin_markers[bin_id].isdisjoint(contig_marker_set)
        ]
        if not possible_bins:
            continue
        work_items.append(
            (
                contig_id,
                possible_bins,
                normalized_tetramer_profiles[contig_id],
                coverages[contig_id],
            )
        )

    if nthreads == 1:
        results = [
            _score_contig_against_bins_exact(
                contig_id=contig_id,
                possible_bins=possible_bins,
                contig_tetramers=contig_tetramers,
                contig_coverage=contig_coverage,
                bin_tetramer_matrices=bin_tetramer_matrices,
                bin_coverage_matrices=bin_coverage_matrices,
                w_intra=w_intra,
            )
            for (
                contig_id,
                possible_bins,
                contig_tetramers,
                contig_coverage,
            ) in work_items
        ]
    elif work_items:
        chunksize = max(1, math.ceil(len(work_items) / (nthreads * 4)))
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=nthreads,
            initializer=_init_contig_scoring_worker,
            initargs=(
                bin_tetramer_matrices,
                bin_coverage_matrices,
                w_intra,
            ),
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
        contig_id, best_bin, _ = result
        # Guard: contig may appear in multiple work items
        if contig_id in bin_of_contig:
            continue
        bins[best_bin].append(contig_id)
        bin_of_contig[contig_id] = best_bin
        binned_contigs_with_markers.append(contig_id)
        bin_markers[best_bin].update(contig_markers[contig_id])

    return bins, bin_of_contig, n_bins, bin_markers, binned_contigs_with_markers
