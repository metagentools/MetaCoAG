#!/usr/bin/env python3

import hashlib
import re

from collections import defaultdict, deque

from Bio import SeqIO

from metacoag.metacoag_utils.bidirectionalmap import BidirectionalMap


def hash_sequence(sequence):
    """Return a stable, compact digest for sequence identity matching."""
    return hashlib.sha256(str(sequence).encode("utf-8")).digest()


def get_segment_paths_spades(contig_paths):
    paths = {}
    segment_contigs = {}
    node_count = 0

    contig_map = BidirectionalMap()
    contig_names = BidirectionalMap()
    current_contig_id = ""

    with open(contig_paths) as file:
        name = file.readline().strip()
        path = file.readline().strip()

        while name != "" and path != "":
            while ";" in path:
                path = path[:-2] + "," + file.readline()

            match = re.search(r"NODE_(.*)_length_", name)
            contig_id = str(int(match.group(1)))

            segments = path.rstrip().split(",")

            if current_contig_id != contig_id:
                contig_map[node_count] = int(contig_id)
                contig_names[node_count] = name.strip()
                current_contig_id = contig_id
                node_count += 1

            if contig_id not in paths:
                paths[contig_id] = segments

            for segment in segments:
                segment_contigs.setdefault(segment, set()).add(contig_id)

            name = file.readline().strip()
            path = file.readline().strip()

    return paths, segment_contigs, node_count, contig_map, contig_names


def get_graph_edges_spades(
    assembly_graph_file, contigs_map, contigs_map_rev, paths, segment_contigs
):
    links = []
    links_map = defaultdict(set)

    # Get links from assembly_graph_with_scaffolds.gfa
    with open(assembly_graph_file) as file:
        line = file.readline()

        while line != "":
            # Identify lines with link information
            if "L" in line:
                fields = line.split("\t")
                first_segment = fields[1] + fields[2]
                second_segment = fields[3] + fields[4]
                links_map[first_segment].add(second_segment)
                links_map[second_segment].add(first_segment)
                links.append(f"{first_segment} {second_segment}")
            line = file.readline()

    # Create list of edges
    edge_list = []

    for contig_id in range(len(paths)):
        segments = paths[str(contigs_map[contig_id])]

        new_links = []

        for segment in segments:
            active_segment = segment
            if active_segment.endswith("+"):
                reverse_segment = active_segment[:-1] + "-"
            else:
                reverse_segment = active_segment[:-1] + "+"

            if segment in links_map:
                new_links.extend(links_map[segment])

            if reverse_segment in links_map:
                new_links.extend(links_map[reverse_segment])

        if active_segment in segment_contigs:
            for linked_contig in segment_contigs[active_segment]:
                linked_contig_id = contigs_map_rev[int(linked_contig)]
                if contig_id != linked_contig_id:
                    # Add edge to list of edges
                    edge_list.append((contig_id, linked_contig_id))

        if reverse_segment in segment_contigs:
            for linked_contig in segment_contigs[reverse_segment]:
                linked_contig_id = contigs_map_rev[int(linked_contig)]
                if contig_id != linked_contig_id:
                    # Add edge to list of edges
                    edge_list.append((contig_id, linked_contig_id))

        for new_link in new_links:
            if new_link in segment_contigs:
                for linked_contig in segment_contigs[new_link]:
                    linked_contig_id = contigs_map_rev[int(linked_contig)]
                    if contig_id != linked_contig_id:
                        # Add edge to list of edges
                        edge_list.append((contig_id, linked_contig_id))

    return edge_list


def get_flye_contig_map(contigs_file):
    contig_names = BidirectionalMap()

    for contig_id, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
        contig_names[contig_id] = record.id

    return contig_names


def get_links_flye(contig_paths, contig_names_rev):
    paths = {}
    segment_contigs = {}

    contig_map = BidirectionalMap()

    with open(contig_paths) as file:
        for line in file:
            if not line.startswith("#"):
                fields = line.strip().split()

                contig_name = fields[0]

                path = fields[-1]
                path = path.replace("*", "")

                if path.startswith(","):
                    path = path[1:]

                if path.endswith(","):
                    path = path[:-1]

                segments = path.rstrip().split(",")

                contig_id = contig_names_rev[contig_name]

                if contig_id not in paths:
                    paths[contig_id] = segments

                for segment in segments:
                    segment_contigs.setdefault(segment, set()).add(contig_id)

    return paths, segment_contigs, len(contig_names_rev), contig_map


def get_graph_edges_flye(
    assembly_graph_file, contigs_map, contigs_map_rev, paths, segment_contigs
):
    links_map = defaultdict(set)

    # Get links from assembly_graph_with_scaffolds.gfa
    with open(assembly_graph_file) as file:
        line = file.readline()

        while line != "":
            # Identify lines with link information
            if "L" in line:
                strings = line.split("\t")

                f1, f2 = "", ""

                if strings[2] == "+":
                    f1 = strings[1][5:]
                if strings[2] == "-":
                    f1 = "-" + strings[1][5:]
                if strings[4] == "+":
                    f2 = strings[3][5:]
                if strings[4] == "-":
                    f2 = "-" + strings[3][5:]

                links_map[f1].add(f2)
                links_map[f2].add(f1)

            line = file.readline()

    # Create list of edges
    edge_list = []

    for i in paths:
        segments = paths[i]

        new_links = []

        for segment in segments:
            my_segment = segment
            my_segment_num = ""

            my_segment_rev = ""

            if my_segment.startswith("-"):
                my_segment_rev = my_segment[1:]
                my_segment_num = my_segment[1:]
            else:
                my_segment_rev = "-" + my_segment
                my_segment_num = my_segment

            if my_segment in links_map:
                new_links.extend(list(links_map[my_segment]))

            if my_segment_rev in links_map:
                new_links.extend(list(links_map[my_segment_rev]))

            if my_segment in segment_contigs:
                for contig in segment_contigs[my_segment]:
                    if i != contig:
                        # Add edge to list of edges
                        edge_list.append((i, contig))

            if my_segment_rev in segment_contigs:
                for contig in segment_contigs[my_segment_rev]:
                    if i != contig:
                        # Add edge to list of edges
                        edge_list.append((i, contig))

            if my_segment_num in segment_contigs:
                for contig in segment_contigs[my_segment_num]:
                    if i != contig:
                        # Add edge to list of edges
                        edge_list.append((i, contig))

        for new_link in new_links:
            if new_link in segment_contigs:
                for contig in segment_contigs[new_link]:
                    if i != contig:
                        # Add edge to list of edges
                        edge_list.append((i, contig))

            if new_link.startswith("-"):
                if new_link[1:] in segment_contigs:
                    for contig in segment_contigs[new_link[1:]]:
                        if i != contig:
                            # Add edge to list of edges
                            edge_list.append((i, contig))

    return edge_list


def get_links_megahit(assembly_graph_file):
    node_count = 0

    graph_contig_hashes = {}

    links = []

    contig_names = BidirectionalMap()

    # Get links from .gfa file
    with open(assembly_graph_file) as file:
        line = file.readline()

        while line != "":
            # Identify lines with link information
            if line.startswith("L"):
                link = []

                fields = line.split("\t")

                link1 = fields[1]
                link2 = fields[3]

                link.append(link1)
                link.append(link2)
                links.append(link)

            elif line.startswith("S"):
                fields = line.split()

                contig_names[node_count] = fields[1]
                graph_contig_hashes[fields[1]] = hash_sequence(fields[2])

                node_count += 1

            line = file.readline()

    return node_count, graph_contig_hashes, links, contig_names


def get_links_megahit_custom(assembly_graph_file):
    contig_names = BidirectionalMap()

    node_count = 0

    links = []

    # Get contig connections from .gfa file
    with open(assembly_graph_file) as file:
        for line in file:
            line = line.strip()

            # Count the number of contigs
            if line.startswith("S"):
                fields = line.split("\t")
                node_name = fields[1][:-2]
                contig_names[node_count] = node_name
                node_count += 1

            # Identify lines with link information
            elif line.startswith("L"):
                link = []
                fields = line.split("\t")

                if fields[1] != fields[3]:
                    link.append(fields[1])
                    link.append(fields[3])
                    links.append(link)

    return node_count, links, contig_names


def get_links_custom(assembly_graph_file):
    contig_names = BidirectionalMap()

    node_count = 0

    links = []

    # Get contig connections from .gfa file
    with open(assembly_graph_file) as file:
        for line in file:
            line = line.strip()

            # Count the number of contigs
            if line.startswith("S"):
                fields = line.split("\t")
                node_name = fields[1]
                contig_names[node_count] = node_name
                node_count += 1

            # Identify lines with link information
            elif line.startswith("L"):
                link = []
                fields = line.split("\t")

                if fields[1] != fields[3]:
                    link.append(fields[1])
                    link.append(fields[3])
                    links.append(link)

    return node_count, links, contig_names


def get_graph_edges_megahit(links, contig_names_rev):
    edge_list = []

    # Iterate links
    for link in links:
        # Remove self loops
        if link[0] != link[1]:
            # Add edge to list of edges
            edge_list.append((contig_names_rev[link[0]], contig_names_rev[link[1]]))

    return edge_list


def get_isolated(node_count, assembly_graph):
    isolated = []

    # Get isolated contigs which have no neighbours
    for i in range(node_count):
        neighbours = assembly_graph.neighbors(i, mode="ALL")

        if len(neighbours) == 0:
            isolated.append(i)

    return isolated


def get_non_isolated(node_count, assembly_graph, binned_contigs, nthreads):
    """Return per-node labelled connected components.

    This preserves the pre-connected-components return shape: one list per graph
    node, with component members only for nodes that are already binned.
    ``nthreads`` is kept for API compatibility.
    """
    del nthreads

    binned_contigs = set(binned_contigs)
    return [
        get_connected_component_contigs(contig_id, assembly_graph, binned_contigs)
        for contig_id in range(node_count)
    ]


def get_connected_component_contigs(contig_id, assembly_graph, binned_contigs):
    if contig_id not in binned_contigs:
        return []

    component = []
    visited = {contig_id}
    queue = deque([contig_id])

    while queue:
        active_contig = queue.popleft()
        component.append(active_contig)

        for neighbour in assembly_graph.neighbors(active_contig, mode="ALL"):
            if neighbour not in visited:
                visited.add(neighbour)
                queue.append(neighbour)

    if not any(contig in binned_contigs for contig in component):
        return []

    return component


# Backward-compatible alias for older imports.
get_connected_components = get_connected_component_contigs
