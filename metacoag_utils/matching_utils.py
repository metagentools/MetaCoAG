#!/usr/bin/env python3

import logging
import math
import operator
import sys

import networkx as nx
from scipy.spatial import distance

# Constants set from MaxBin 2.0
MU_INTRA, SIGMA_INTRA = 0, 0.01037897 / 2
MU_INTER, SIGMA_INTER = 0.0676654, 0.03419337
VERY_SMALL_DOUBLE = 1e-10
MAX_WEIGHT = sys.float_info.max

# create logger
logger = logging.getLogger("MetaCoaAG 1.1.2")


def normpdf(x, mean, sd):
    var = float(sd) ** 2
    denom = sd * (2 * math.pi) ** 0.5
    num = math.exp(-((float(x) - float(mean)) ** 2) / (2 * var))
    return num / denom


def get_tetramer_distance(seq1, seq2):
    return distance.euclidean(seq1, seq2)


def get_coverage_distance(cov1, cov2):
    return distance.euclidean(cov1, cov2)


def get_comp_probability(tetramer_dist):
    gaus_intra = normpdf(tetramer_dist, MU_INTRA, SIGMA_INTRA)
    gaus_inter = normpdf(tetramer_dist, MU_INTER, SIGMA_INTER)
    return gaus_intra / (gaus_intra + gaus_inter)


def get_cov_probability(cov1, cov2):
    poisson_prod_1 = 1
    poisson_prod_2 = 1

    for i in range(len(cov1)):
        # Adapted from http://www.masaers.com/2013/10/08/Implementing-Poisson-pmf.html
        poisson_pmf_1 = math.exp(
            (cov1[i] * math.log(cov2[i])) - math.lgamma(cov1[i] + 1.0) - cov2[i]
        )

        poisson_pmf_2 = math.exp(
            (cov2[i] * math.log(cov1[i])) - math.lgamma(cov2[i] + 1.0) - cov1[i]
        )

        if poisson_pmf_1 < VERY_SMALL_DOUBLE:
            poisson_pmf_1 = VERY_SMALL_DOUBLE

        if poisson_pmf_2 < VERY_SMALL_DOUBLE:
            poisson_pmf_2 = VERY_SMALL_DOUBLE

        poisson_prod_1 = poisson_prod_1 * poisson_pmf_1

        poisson_prod_2 = poisson_prod_2 * poisson_pmf_2

    return min(poisson_prod_1, poisson_prod_2)


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

            if len(to_bin) != 0:
                for contig in to_bin:
                    contigid = contig

                    if contigid not in top_nodes:
                        top_nodes.append(contigid)

                    for b in range(n_bins):
                        log_prob_sum = 0
                        n_contigs = len(bins[b])

                        for j in range(n_contigs):
                            tetramer_dist = get_tetramer_distance(
                                normalized_tetramer_profiles[contigid],
                                normalized_tetramer_profiles[bins[b][j]],
                            )
                            prob_comp = get_comp_probability(tetramer_dist)
                            prob_cov = get_cov_probability(
                                coverages[contigid], coverages[bins[b][j]]
                            )

                            prob_product = prob_comp * prob_cov

                            log_prob = 0

                            if prob_product > 0.0:
                                log_prob = -(
                                    math.log(prob_comp, 10) + math.log(prob_cov, 10)
                                )
                            else:
                                log_prob = MAX_WEIGHT

                            log_prob_sum += log_prob

                        if log_prob_sum != float("inf"):
                            edges.append(
                                (bins[b][0], contigid, log_prob_sum / n_contigs)
                            )
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
                                path_len_sum = 0

                                for contig_in_bin in bins[b]:
                                    shortest_paths = assembly_graph.get_shortest_paths(
                                        my_matching[l], to=contig_in_bin
                                    )

                                    if len(shortest_paths) != 0:
                                        path_len_sum += len(shortest_paths[0])

                                avg_path_len = math.floor(path_len_sum / len(bins[b]))

                                # logger.debug("To assign contig " + contig_names[my_matching[l]] + " to bin "+str(
                                #     b+1) + " based on contig " + str(l) + " weight="+str(edge_weights[(l, my_matching[l])]))

                                if (
                                    edge_weights[(l, my_matching[l])] <= w_intra
                                    and avg_path_len <= d_limit
                                ):
                                    can_assign = False

                                    common_mgs = set(bin_markers[b]).intersection(
                                        set(contig_markers[my_matching[l]])
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

                                        bin_markers[b] = list(
                                            set(
                                                bin_markers[b]
                                                + contig_markers[my_matching[l]]
                                            )
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
                        path_len_sum = 0

                        for contig_in_bin in bins[not_binned[longest_nb_contig][1]]:
                            shortest_paths = assembly_graph.get_shortest_paths(
                                longest_nb_contig, to=contig_in_bin
                            )

                            if len(shortest_paths) != 0:
                                path_len_sum += len(shortest_paths[0])

                        avg_path_len = path_len_sum / len(
                            bins[not_binned[longest_nb_contig][1]]
                        )

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

                            bin_markers[n_bins] = contig_markers[longest_nb_contig]
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
):
    for contig in unbinned_mg_contigs:
        if contig[1] >= min_length:
            possible_bins = []

            for b in bin_markers:
                common_mgs = list(
                    set(bin_markers[b]).intersection(set(contig_markers[contig[0]]))
                )
                if len(common_mgs) == 0:
                    possible_bins.append(b)

            if len(possible_bins) != 0:
                contigid = contig[0]

                bin_weights = []

                for b in possible_bins:
                    log_prob_sum = 0
                    n_contigs = len(bins[b])

                    for j in range(n_contigs):
                        tetramer_dist = get_tetramer_distance(
                            normalized_tetramer_profiles[contigid],
                            normalized_tetramer_profiles[bins[b][j]],
                        )
                        prob_comp = get_comp_probability(tetramer_dist)
                        prob_cov = get_cov_probability(
                            coverages[contigid], coverages[bins[b][j]]
                        )

                        prob_product = prob_comp * prob_cov

                        log_prob = 0

                        if prob_product > 0.0:
                            log_prob = -(
                                math.log(prob_comp, 10) + math.log(prob_cov, 10)
                            )
                        else:
                            log_prob = MAX_WEIGHT

                        log_prob_sum += log_prob

                    if log_prob_sum != float("inf"):
                        bin_weights.append(log_prob_sum / n_contigs)
                    else:
                        bin_weights.append(MAX_WEIGHT)

                min_b_index, min_b_value = min(
                    enumerate(bin_weights), key=operator.itemgetter(1)
                )

                if min_b_value <= w_intra:
                    bins[possible_bins[min_b_index]].append(contigid)
                    bin_of_contig[contigid] = possible_bins[min_b_index]
                    binned_contigs_with_markers.append(contigid)

                    bin_markers[possible_bins[min_b_index]] = list(
                        set(
                            bin_markers[possible_bins[min_b_index]]
                            + contig_markers[contigid]
                        )
                    )

    return bins, bin_of_contig, n_bins, bin_markers, binned_contigs_with_markers
