#!/usr/bin/env python3

import concurrent.futures
import csv
import logging
import math
import operator
import pathlib
import shutil
import sys

from collections import OrderedDict
from dataclasses import dataclass
from typing import Any, Optional

from Bio import SeqIO
from igraph import Graph
from tqdm import tqdm

from metacoag.metacoag_utils import (
    feature_utils,
    graph_utils,
    label_prop_utils,
    marker_gene_utils,
    matching_utils,
)
from metacoag.metacoag_utils.bidirectionalmap import BidirectionalMap


MAX_WEIGHT = sys.float_info.max
MAX_OPEN_BIN_FILES = 32


@dataclass(frozen=True)
class PipelineConfig:
    assembler: str
    graph: str
    contigs: str
    abundance: str
    paths: Optional[str]
    output: pathlib.Path
    hmm: str
    prefix: str
    min_length: int
    p_intra: float
    p_inter: float
    d_limit: int
    depth: int
    n_mg: int
    no_cut_tc: bool
    mg_threshold: float
    bin_mg_threshold: float
    min_bin_size: int
    delimiter: str
    nthreads: int

    @classmethod
    def from_args(cls, args):
        prefix = args.prefix
        if prefix and not prefix.endswith("_"):
            prefix = f"{prefix}_"

        return cls(
            assembler=args.assembler.lower(),
            graph=args.graph,
            contigs=args.contigs,
            abundance=args.abundance,
            paths=args.paths,
            output=pathlib.Path(args.output),
            hmm=args.hmm,
            prefix=prefix,
            min_length=args.min_length,
            p_intra=args.p_intra,
            p_inter=args.p_inter,
            d_limit=args.d_limit,
            depth=args.depth,
            n_mg=args.n_mg,
            no_cut_tc=args.no_cut_tc,
            mg_threshold=args.mg_threshold,
            bin_mg_threshold=args.bin_mg_threshold,
            min_bin_size=args.min_bin_size,
            delimiter=args.delimiter,
            nthreads=args.nthreads,
        )


@dataclass
class GraphData:
    node_count: int
    contig_names: Any
    contig_names_rev: Any
    assembly_graph: Graph
    graph_to_contig_map: Any = None
    graph_to_contig_map_rev: Any = None
    contig_descriptions: Optional[dict] = None


@dataclass
class FeatureData:
    coverages: Any
    contig_lengths: Any
    normalized_tetramer_profiles: Any
    n_samples: int
    w_intra: float
    w_inter: float


@dataclass
class MarkerData:
    marker_contigs: dict
    marker_contig_counts: dict
    contig_markers: dict


@dataclass
class BinState:
    bins: dict
    bin_of_contig: dict
    bin_markers: dict
    binned_contigs_with_markers: list
    n_bins: int
    smg_bin_counts: list
    seed_tetramer_profiles: dict
    seed_coverage_profiles: dict


@dataclass
class MergePlan:
    bin_cliques: list
    bin_clique_sizes: dict
    bins_to_remove: set


class BinFastaWriter:
    def __init__(self, paths, max_open_files=MAX_OPEN_BIN_FILES):
        if max_open_files <= 0:
            raise ValueError("max_open_files must be positive")

        self.paths = set(paths)
        self.max_open_files = max_open_files
        self.open_files = OrderedDict()

    def __enter__(self):
        for path in self.paths:
            with open(path, "w"):
                pass
        return self

    def write(self, path, record):
        output_file = self.open_files.pop(path, None)
        if output_file is None:
            if len(self.open_files) >= self.max_open_files:
                _, oldest_file = self.open_files.popitem(last=False)
                oldest_file.close()
            output_file = open(path, "a")

        self.open_files[path] = output_file
        output_file.write(f">{record.id}\n{record.seq}\n")

    def __exit__(self, exc_type, exc_value, traceback):
        for output_file in self.open_files.values():
            output_file.close()
        self.open_files.clear()


def validate_inputs(config):
    required_files = [config.graph, config.contigs, config.abundance]
    if config.assembler in {"spades", "flye"}:
        if config.paths is None:
            expected = (
                "contigs.paths" if config.assembler == "spades" else "assembly_info.txt"
            )
            raise ValueError(f"{config.assembler} requires a {expected} file")
        required_files.append(config.paths)

    missing_files = [
        path for path in required_files if not pathlib.Path(path).is_file()
    ]
    if missing_files:
        raise FileNotFoundError(f"Input file does not exist: {missing_files[0]}")

    positive_values = {
        "min_bin_size": config.min_bin_size,
        "depth": config.depth,
        "d_limit": config.d_limit,
        "nthreads": config.nthreads,
    }
    for name, value in positive_values.items():
        if value <= 0:
            raise ValueError(f"{name} must be positive")

    if not 0 < config.p_intra <= 1:
        raise ValueError("p_intra must be greater than 0 and at most 1")
    if not 0 < config.p_inter <= 1:
        raise ValueError("p_inter must be greater than 0 and at most 1")

    config.output.mkdir(parents=True, exist_ok=True)


def log_inputs(config, logger):
    logger.info(
        "Welcome to MetaCoAG: Binning Metagenomic Contigs via Composition, "
        "Coverage and Assembly Graphs."
    )
    logger.info("Input arguments:")
    logger.info("Assembler used: %s", config.assembler)
    logger.info("Contigs file: %s", config.contigs)
    logger.info("Assembly graph file: %s", config.graph)
    logger.info("Contig paths file: %s", config.paths)
    logger.info("Abundance file: %s", config.abundance)
    logger.info("Final binning output file: %s", config.output)
    logger.info("Marker gene file hmm: %s", config.hmm)
    logger.debug("Number of marker genes in hmm file: %s", config.n_mg)
    logger.info("Minimum length of contigs to consider: %s", config.min_length)
    logger.info("Depth to consider for label propagation: %s", config.depth)
    logger.info("p_intra: %s", config.p_intra)
    logger.info("p_inter: %s", config.p_inter)
    logger.info("Do not use --cut_tc: %s", config.no_cut_tc)
    logger.info("mg_threshold: %s", config.mg_threshold)
    logger.info("bin_mg_threshold: %s", config.bin_mg_threshold)
    logger.info("min_bin_size: %s base pairs", config.min_bin_size)
    logger.info("d_limit: %s", config.d_limit)
    logger.info("Number of threads: %s", config.nthreads)


def load_graph(config, logger):
    assembler = config.assembler
    graph_to_contig_map = None
    graph_to_contig_map_rev = None
    contig_descriptions = None

    if assembler == "spades":
        (
            paths,
            segment_contigs,
            node_count,
            contigs_map,
            contig_names,
        ) = graph_utils.get_segment_paths_spades(config.paths)
        contigs_map_rev = contigs_map.inverse
        contig_names_rev = contig_names.inverse
    elif assembler == "megahit":
        original_hash_to_name = {}
        contig_descriptions = {}
        for record in SeqIO.parse(config.contigs, "fasta"):
            original_hash_to_name[graph_utils.hash_sequence(record.seq)] = record.id
            contig_descriptions[record.id] = record.description

        (
            node_count,
            graph_contig_hashes,
            links,
            contig_names,
        ) = graph_utils.get_links_megahit(config.graph)
        contig_names_rev = contig_names.inverse
    elif assembler == "flye":
        contig_names = graph_utils.get_flye_contig_map(config.contigs)
        contig_names_rev = contig_names.inverse
        (
            paths,
            segment_contigs,
            node_count,
            contigs_map,
        ) = graph_utils.get_links_flye(config.paths, contig_names_rev)
        contigs_map_rev = contigs_map.inverse
    elif assembler == "megahitc":
        node_count, links, contig_names = graph_utils.get_links_megahit_custom(
            config.graph
        )
        contig_names_rev = contig_names.inverse
    elif assembler == "custom":
        node_count, links, contig_names = graph_utils.get_links_custom(config.graph)
        contig_names_rev = contig_names.inverse
    else:
        raise ValueError(f"Unsupported assembler: {assembler}")

    assembly_graph = Graph()
    assembly_graph.add_vertices(node_count)
    logger.info("Total number of contigs available: %s", node_count)

    for i in range(node_count):
        assembly_graph.vs[i]["id"] = i
        assembly_graph.vs[i]["label"] = contig_names[i]

    if assembler == "spades":
        edge_list = graph_utils.get_graph_edges_spades(
            assembly_graph_file=config.graph,
            contigs_map=contigs_map,
            contigs_map_rev=contigs_map_rev,
            paths=paths,
            segment_contigs=segment_contigs,
        )
    elif assembler == "flye":
        edge_list = graph_utils.get_graph_edges_flye(
            assembly_graph_file=config.graph,
            contigs_map=contigs_map,
            contigs_map_rev=contigs_map_rev,
            paths=paths,
            segment_contigs=segment_contigs,
        )
    else:
        edge_list = graph_utils.get_graph_edges_megahit(
            links=links, contig_names_rev=contig_names_rev
        )

    assembly_graph.add_edges(edge_list)
    assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)
    logger.info(
        "Total number of edges in the assembly graph: %s", assembly_graph.ecount()
    )

    if assembler == "megahit":
        graph_to_contig_map = BidirectionalMap()
        for graph_name, graph_hash in graph_contig_hashes.items():
            if graph_hash in original_hash_to_name:
                graph_to_contig_map[graph_name] = original_hash_to_name[graph_hash]
        graph_to_contig_map_rev = graph_to_contig_map.inverse

    isolated = graph_utils.get_isolated(node_count, assembly_graph)
    logger.info("Total isolated contigs in the assembly graph: %s", len(isolated))

    graph_data = GraphData(
        node_count=node_count,
        contig_names=contig_names,
        contig_names_rev=contig_names_rev,
        assembly_graph=assembly_graph,
        graph_to_contig_map=graph_to_contig_map,
        graph_to_contig_map_rev=graph_to_contig_map_rev,
        contig_descriptions=contig_descriptions,
    )
    return graph_data, isolated


def extract_features(config, graph_data, isolated, logger):
    logger.info("Obtaining lengths and coverage values of contigs")

    if config.assembler == "megahit":
        coverages, contig_lengths, n_samples = feature_utils.get_cov_len_megahit(
            contigs_file=config.contigs,
            contig_names_rev=graph_data.contig_names_rev,
            graph_to_contig_map_rev=graph_data.graph_to_contig_map_rev,
            min_length=config.min_length,
            abundance_file=config.abundance,
        )
    else:
        coverages, contig_lengths, n_samples = feature_utils.get_cov_len(
            contigs_file=config.contigs,
            contig_names_rev=graph_data.contig_names_rev,
            min_length=config.min_length,
            abundance_file=config.abundance,
        )

    long_contig_count = int((contig_lengths >= config.min_length).sum())
    isolated_long_count = sum(
        contig_lengths[contig] >= config.min_length for contig in isolated
    )
    logger.info("Total long contigs: %s", long_contig_count)
    logger.info(
        "Total isolated long contigs in the assembly graph: %s",
        isolated_long_count,
    )

    bin_threshold = -math.log(config.p_intra, 10)
    break_threshold = -math.log(config.p_inter, 10)
    w_intra = bin_threshold * (n_samples + 1)
    w_inter = break_threshold * (n_samples + 1)
    logger.debug("w_intra: %s", w_intra)
    logger.debug("w_inter: %s", w_inter)

    logger.info("Obtaining tetranucleotide frequencies of contigs")
    normalized_tetramer_profiles = feature_utils.get_tetramer_profiles(
        output_path=config.output,
        contigs_file=config.contigs,
        contig_names_rev=graph_data.contig_names_rev,
        contig_lengths=contig_lengths,
        min_length=config.min_length,
        nthreads=config.nthreads,
        graph_to_contig_map_rev=graph_data.graph_to_contig_map_rev,
    )

    return FeatureData(
        coverages=coverages,
        contig_lengths=contig_lengths,
        normalized_tetramer_profiles=normalized_tetramer_profiles,
        n_samples=n_samples,
        w_intra=w_intra,
        w_inter=w_inter,
    )


def parse_markers(config, graph_data, features, logger):
    logger.info("Scanning for single-copy marker genes")
    hmm_output = pathlib.Path(f"{config.contigs}.hmmout")

    if not hmm_output.exists():
        missing_tools = [
            tool
            for tool in ("run_FragGeneScan.pl", "hmmsearch")
            if shutil.which(tool) is None
        ]
        if missing_tools:
            raise FileNotFoundError(
                "Required marker-gene tool is not installed: "
                + ", ".join(missing_tools)
            )

        logger.info("Obtaining hmmout file")
        marker_gene_utils.scan_for_marker_genes(
            contigs_file=config.contigs,
            nthreads=config.nthreads,
            marker_url=config.hmm,
            no_cut_tc=config.no_cut_tc,
        )
    else:
        logger.info(".hmmout file already exists")

    logger.info("Obtaining contigs with single-copy marker genes")
    if config.assembler == "megahit":
        marker_values = marker_gene_utils.get_contigs_with_marker_genes_megahit(
            contigs_file=config.contigs,
            contig_names_rev=graph_data.contig_names_rev,
            graph_to_contig_map_rev=graph_data.graph_to_contig_map_rev,
            mg_length_threshold=config.mg_threshold,
            contig_lengths=features.contig_lengths,
            min_length=config.min_length,
        )
    else:
        marker_values = marker_gene_utils.get_contigs_with_marker_genes(
            contigs_file=config.contigs,
            contig_names_rev=graph_data.contig_names_rev,
            mg_length_threshold=config.mg_threshold,
            contig_lengths=features.contig_lengths,
            min_length=config.min_length,
        )

    marker_data = MarkerData(*marker_values)
    logger.info(
        "Number of contigs containing single-copy marker genes: %s",
        len(marker_data.contig_markers),
    )
    if not marker_data.contig_markers:
        raise RuntimeError(
            "No contigs contain qualifying single-copy marker genes; "
            "the dataset cannot be binned"
        )

    return marker_data


def _build_marker_iterations(marker_data):
    marker_counts = sorted(marker_data.marker_contig_counts.values(), reverse=True)
    iterations = {}
    iteration_index = 0

    for count in sorted(set(marker_counts), reverse=True):
        marker_scores = {}
        for marker, marker_count in marker_data.marker_contig_counts.items():
            if marker_count == count:
                marker_scores[marker] = sum(
                    len(marker_data.contig_markers[contig])
                    for contig in marker_data.marker_contigs[marker]
                )

        for marker, _ in sorted(
            marker_scores.items(), key=operator.itemgetter(1), reverse=True
        ):
            iterations[iteration_index] = marker_data.marker_contigs[marker]
            iteration_index += 1

    return iterations, marker_counts


def match_seed_bins(config, graph_data, features, marker_data, logger):
    logger.info("Determining contig counts for each single-copy marker gene")
    smg_iteration, marker_counts = _build_marker_iterations(marker_data)
    logger.debug("Contig counts of single-copy marker genes: %s", marker_counts)

    bins = {}
    bin_of_contig = {}
    bin_markers = {}
    binned_contigs_with_markers = []

    logger.info("Initialising bins")
    for bin_id, contig_num in enumerate(smg_iteration[0]):
        binned_contigs_with_markers.append(contig_num)
        bins[bin_id] = [contig_num]
        bin_of_contig[contig_num] = bin_id
        bin_markers[bin_id] = marker_data.contig_markers[contig_num].copy()

    logger.debug("Number of initial bins detected: %s", len(smg_iteration[0]))
    logger.info("Matching and assigning contigs with single-copy marker genes to bins")
    (
        bins,
        bin_of_contig,
        n_bins,
        bin_markers,
        binned_contigs_with_markers,
    ) = matching_utils.match_contigs(
        smg_iteration=smg_iteration,
        bins=bins,
        n_bins=0,
        bin_of_contig=bin_of_contig,
        binned_contigs_with_markers=binned_contigs_with_markers,
        bin_markers=bin_markers,
        contig_markers=marker_data.contig_markers,
        contig_lengths=features.contig_lengths,
        contig_names=graph_data.contig_names,
        normalized_tetramer_profiles=features.normalized_tetramer_profiles,
        coverages=features.coverages,
        assembly_graph=graph_data.assembly_graph,
        w_intra=features.w_intra,
        w_inter=features.w_inter,
        d_limit=config.d_limit,
    )
    logger.debug("Number of bins after matching: %s", len(bins))

    unbinned_marker_contigs = set(marker_data.contig_markers) - set(
        binned_contigs_with_markers
    )
    unbinned_by_length = sorted(
        (
            (contig, features.contig_lengths[contig])
            for contig in unbinned_marker_contigs
        ),
        key=operator.itemgetter(1),
        reverse=True,
    )
    logger.debug(
        "Number of unbinned contigs with single-copy marker genes: %s",
        len(unbinned_marker_contigs),
    )
    logger.info("Further assigning contigs with single-copy marker genes")

    (
        bins,
        bin_of_contig,
        n_bins,
        bin_markers,
        binned_contigs_with_markers,
    ) = matching_utils.further_match_contigs(
        unbinned_mg_contigs=unbinned_by_length,
        min_length=config.min_length,
        bins=bins,
        n_bins=n_bins,
        bin_of_contig=bin_of_contig,
        binned_contigs_with_markers=binned_contigs_with_markers,
        bin_markers=bin_markers,
        contig_markers=marker_data.contig_markers,
        normalized_tetramer_profiles=features.normalized_tetramer_profiles,
        coverages=features.coverages,
        w_intra=features.w_intra,
        nthreads=config.nthreads,
    )

    remaining_marker_contigs = set(marker_data.contig_markers) - set(
        binned_contigs_with_markers
    )
    logger.debug(
        "Remaining number of unbinned MG seed contigs: %s",
        len(remaining_marker_contigs),
    )
    logger.debug(
        "Number of binned contigs with single-copy marker genes: %s",
        len(bin_of_contig),
    )

    smg_bin_counts = [len(bins[bin_id]) for bin_id in bins]
    seed_tetramer_profiles, seed_coverage_profiles = feature_utils.get_bin_profiles(
        bins=bins,
        coverages=features.coverages,
        normalized_tetramer_profiles=features.normalized_tetramer_profiles,
    )

    return BinState(
        bins=bins,
        bin_of_contig=bin_of_contig,
        bin_markers=bin_markers,
        binned_contigs_with_markers=binned_contigs_with_markers,
        n_bins=n_bins,
        smg_bin_counts=smg_bin_counts,
        seed_tetramer_profiles=seed_tetramer_profiles,
        seed_coverage_profiles=seed_coverage_profiles,
    )


def _update_bin_state(state, values):
    (
        state.bins,
        state.bin_of_contig,
        state.bin_markers,
        state.binned_contigs_with_markers,
    ) = values


def propagate_bins(config, graph_data, features, marker_data, state, logger):
    binned_contigs = list(state.bin_of_contig)
    non_isolated = graph_utils.get_non_isolated(
        node_count=graph_data.node_count,
        assembly_graph=graph_data.assembly_graph,
        binned_contigs=binned_contigs,
        nthreads=config.nthreads,
    )
    logger.debug("Number of non-isolated contigs: %s", len(non_isolated))
    logger.info("Propagating labels to connected vertices of unlabelled long contigs")

    common_args = {
        "contig_markers": marker_data.contig_markers,
        "smg_bin_counts": state.smg_bin_counts,
        "non_isolated": non_isolated,
        "contig_lengths": features.contig_lengths,
        "min_length": config.min_length,
        "assembly_graph": graph_data.assembly_graph,
        "normalized_tetramer_profiles": features.normalized_tetramer_profiles,
        "coverages": features.coverages,
        "nthreads": config.nthreads,
    }

    _update_bin_state(
        state,
        label_prop_utils.label_prop(
            bin_of_contig=state.bin_of_contig,
            bins=state.bins,
            bin_markers=state.bin_markers,
            binned_contigs_with_markers=state.binned_contigs_with_markers,
            depth=1,
            weight=features.w_intra,
            **common_args,
        ),
    )
    logger.debug("Total number of binned contigs: %s", len(state.bin_of_contig))

    _update_bin_state(
        state,
        label_prop_utils.label_prop(
            bin_of_contig=state.bin_of_contig,
            bins=state.bins,
            bin_markers=state.bin_markers,
            binned_contigs_with_markers=state.binned_contigs_with_markers,
            depth=config.depth,
            weight=features.w_inter,
            **common_args,
        ),
    )
    logger.debug("Total number of binned contigs: %s", len(state.bin_of_contig))

    logger.info("Further propagating labels to vertices of unlabelled long contigs")
    long_unbinned = [
        contig
        for contig in range(graph_data.node_count)
        if contig not in state.bin_of_contig
        and features.contig_lengths[contig] >= config.min_length
    ]

    def assign_long(contig):
        return label_prop_utils.assign_long(
            contig_id=contig,
            coverages=features.coverages,
            normalized_tetramer_profiles=features.normalized_tetramer_profiles,
            bin_tetramer_profiles=state.seed_tetramer_profiles,
            bin_coverage_profiles=state.seed_coverage_profiles,
        )

    with concurrent.futures.ThreadPoolExecutor(max_workers=config.nthreads) as executor:
        assigned = list(
            tqdm(
                executor.map(assign_long, long_unbinned),
                total=len(long_unbinned),
            )
        )

    put_to_bins = [assignment for assignment in assigned if assignment is not None]
    if put_to_bins:
        _update_bin_state(
            state,
            label_prop_utils.assign_to_bins(
                put_to_bins=put_to_bins,
                bins=state.bins,
                bin_of_contig=state.bin_of_contig,
                bin_markers=state.bin_markers,
                binned_contigs_with_markers=state.binned_contigs_with_markers,
                contig_markers=marker_data.contig_markers,
                contig_lengths=features.contig_lengths,
            ),
        )
    else:
        logger.debug("No further contigs were binned")

    logger.info(
        "Further propagating labels to connected vertices of unlabelled long contigs"
    )
    _update_bin_state(
        state,
        label_prop_utils.final_label_prop(
            bin_of_contig=state.bin_of_contig,
            bins=state.bins,
            contig_markers=marker_data.contig_markers,
            bin_markers=state.bin_markers,
            binned_contigs_with_markers=state.binned_contigs_with_markers,
            smg_bin_counts=state.smg_bin_counts,
            contig_lengths=features.contig_lengths,
            min_length=config.min_length,
            assembly_graph=graph_data.assembly_graph,
            normalized_tetramer_profiles=features.normalized_tetramer_profiles,
            coverages=features.coverages,
            depth=config.depth,
            weight=MAX_WEIGHT,
            nthreads=config.nthreads,
        ),
    )
    logger.debug("Total number of binned contigs: %s", len(state.bin_of_contig))
    return state


def merge_bins(config, features, state, logger):
    bin_sizes = {
        bin_id: sum(features.contig_lengths[contig] for contig in members)
        for bin_id, members in state.bins.items()
    }

    bins_graph = Graph()
    bins_graph.add_vertices(len(state.bins))
    for i in range(len(state.bins)):
        bins_graph.vs[i]["id"] = i
        bins_graph.vs[i]["label"] = f"bin {i + 1}"

    bin_ids = list(state.bins)
    marker_sets = {bin_id: frozenset(state.bin_markers[bin_id]) for bin_id in bin_ids}
    best_bins = {bin_id: -1 for bin_id in bin_ids}
    best_weights = {bin_id: MAX_WEIGHT for bin_id in bin_ids}

    for bin_id in bin_ids:
        logger.debug(
            "Bin %s: # contigs: %s, bin size: %sbp, # markers: %s",
            bin_id,
            len(state.bins[bin_id]),
            bin_sizes[bin_id],
            len(marker_sets[bin_id]),
        )

    for index, bin_id in enumerate(bin_ids):
        for other_index in range(index + 1, len(bin_ids)):
            other_bin = bin_ids[other_index]
            if not marker_sets[bin_id].isdisjoint(marker_sets[other_bin]):
                continue

            tetramer_dist = matching_utils.get_tetramer_distance(
                state.seed_tetramer_profiles[bin_id],
                state.seed_tetramer_profiles[other_bin],
            )
            prob_comp = matching_utils.get_comp_probability(tetramer_dist)
            prob_cov = matching_utils.get_cov_probability(
                state.seed_coverage_profiles[bin_id],
                state.seed_coverage_profiles[other_bin],
            )
            if prob_comp * prob_cov > 0.0:
                weight = -(math.log(prob_comp, 10) + math.log(prob_cov, 10))
            else:
                weight = MAX_WEIGHT

            if weight > features.w_intra:
                continue

            if weight < best_weights[bin_id]:
                best_weights[bin_id] = weight
                best_bins[bin_id] = other_bin
            if weight < best_weights[other_bin]:
                best_weights[other_bin] = weight
                best_bins[other_bin] = bin_id

    bins_to_remove = set()
    for bin_id in bin_ids:
        if best_bins[bin_id] != -1:
            bins_graph.add_edge(bin_id, best_bins[bin_id])
        elif len(marker_sets[bin_id]) < config.n_mg * config.bin_mg_threshold:
            bins_to_remove.add(bin_id)

    bin_cliques = bins_graph.maximal_cliques()
    clique_sizes = {}
    for clique in bin_cliques:
        clique_name = "_".join(str(bin_id) for bin_id in clique)
        clique_sizes[clique_name] = sum(bin_sizes[bin_id] for bin_id in clique)

    return MergePlan(
        bin_cliques=bin_cliques,
        bin_clique_sizes=clique_sizes,
        bins_to_remove=bins_to_remove,
    )


def write_output(config, graph_data, state, merge_plan, logger):
    output_bins_path = config.output / f"{config.prefix}bins"
    low_quality_path = config.output / f"{config.prefix}low_quality_bins"
    output_bins_path.mkdir(parents=True, exist_ok=True)
    low_quality_path.mkdir(parents=True, exist_ok=True)

    final_bins = {}
    low_quality_bins = {}
    final_bin_count = 0
    mapping_path = config.output / f"{config.prefix}contig_to_bin.tsv"

    with open(mapping_path, mode="w") as out_file:
        output_writer = csv.writer(
            out_file,
            delimiter=config.delimiter,
            quotechar='"',
            quoting=csv.QUOTE_MINIMAL,
        )

        for clique in merge_plan.bin_cliques:
            bin_name = "_".join(str(bin_id) for bin_id in clique)
            is_removed_singleton = (
                len(clique) == 1 and clique[0] in merge_plan.bins_to_remove
            )
            is_final = (
                not is_removed_singleton
                and merge_plan.bin_clique_sizes[bin_name] >= config.min_bin_size
            )

            if is_final:
                final_bin_count += 1
                destination = final_bins
            else:
                destination = low_quality_bins

            for bin_id in clique:
                for contig in state.bins[bin_id]:
                    destination[contig] = bin_name
                    if not is_final:
                        continue

                    if config.assembler == "megahit":
                        original_name = graph_data.graph_to_contig_map[
                            graph_data.contig_names[contig]
                        ]
                        output_writer.writerow(
                            [
                                graph_data.contig_descriptions[original_name],
                                f"bin_{bin_name}",
                            ]
                        )
                    else:
                        output_writer.writerow(
                            [graph_data.contig_names[contig], f"bin_{bin_name}"]
                        )

    logger.info("Writing the Final Binning result to file")
    final_paths = {
        bin_name: output_bins_path / f"{config.prefix}bin_{bin_name}.fasta"
        for bin_name in set(final_bins.values())
    }
    low_quality_paths = {
        bin_name: low_quality_path / f"{config.prefix}bin_{bin_name}_seqs.fasta"
        for bin_name in set(low_quality_bins.values())
    }

    with BinFastaWriter(
        [*final_paths.values(), *low_quality_paths.values()]
    ) as bin_writer:
        for record in tqdm(
            SeqIO.parse(config.contigs, "fasta"),
            desc="Splitting contigs into bins",
        ):
            if config.assembler == "megahit":
                graph_name = graph_data.graph_to_contig_map_rev[record.id]
                contig_num = graph_data.contig_names_rev[graph_name]
            else:
                contig_num = graph_data.contig_names_rev[record.id]

            if contig_num in final_bins:
                bin_writer.write(final_paths[final_bins[contig_num]], record)
            elif contig_num in low_quality_bins:
                bin_writer.write(
                    low_quality_paths[low_quality_bins[contig_num]], record
                )

    logger.info("Producing %s bins...", final_bin_count)
    logger.info("Final binning results can be found in %s", output_bins_path)
    return final_bin_count
