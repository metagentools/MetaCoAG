#!/usr/bin/env python3

import re

from collections import defaultdict
from metacoag_utils.bidirectionalmap import BidirectionalMap


def get_segment_paths_spades(contig_paths):

    paths = {}
    segment_contigs = {}
    node_count = 0

    my_map = BidirectionalMap()

    contig_names = BidirectionalMap()

    current_contig_num = ""

    contig_segments = {}

    with open(contig_paths) as file:
        name = file.readline().strip()
        path = file.readline().strip()

        while name != "" and path != "":

            while ";" in path:
                path = path[:-2] + "," + file.readline()

            start = 'NODE_'
            end = '_length_'
            contig_num = str(
                int(re.search('%s(.*)%s' % (start, end), name).group(1)))

            segments = path.rstrip().split(",")

            if current_contig_num != contig_num:
                my_map[node_count] = int(contig_num)
                contig_names[node_count] = name.strip()
                current_contig_num = contig_num
                node_count += 1

            if contig_num not in paths:
                paths[contig_num] = [segments[0], segments[-1]]

            for segment in segments:

                if name not in contig_segments:
                    contig_segments[name] = set([segment[:-1]])
                else:
                    contig_segments[name].add(segment[:-1])

                if segment not in segment_contigs:
                    segment_contigs[segment] = set([contig_num])
                else:
                    segment_contigs[segment].add(contig_num)

            name = file.readline().strip()
            path = file.readline().strip()

    return paths, segment_contigs, contig_segments, node_count, my_map, contig_names


def get_graph_edges_spades(
        assembly_graph_file, contigs_map,
        contigs_map_rev, paths, segment_contigs):

    links = []
    links_map = defaultdict(set)

    # Get links from assembly_graph_with_scaffolds.gfa
    with open(assembly_graph_file) as file:
        line = file.readline()

        while line != "":

            # Identify lines with link information
            if "L" in line:
                strings = line.split("\t")
                f1, f2 = strings[1] + strings[2], strings[3] + strings[4]
                links_map[f1].add(f2)
                links_map[f2].add(f1)
                links.append(strings[1] + strings[2] +
                             " " + strings[3] + strings[4])
            line = file.readline()

    # Create list of edges
    edge_list = []

    for i in range(len(paths)):
        segments = paths[str(contigs_map[i])]

        start = segments[0]
        start_rev = ""

        if start.endswith("+"):
            start_rev = start[:-1] + "-"
        else:
            start_rev = start[:-1] + "+"

        end = segments[1]
        end_rev = ""

        if end.endswith("+"):
            end_rev = end[:-1] + "-"
        else:
            end_rev = end[:-1] + "+"

        new_links = []

        if start in links_map:
            new_links.extend(list(links_map[start]))
        if start_rev in links_map:
            new_links.extend(list(links_map[start_rev]))
        if end in links_map:
            new_links.extend(list(links_map[end]))
        if end_rev in links_map:
            new_links.extend(list(links_map[end_rev]))

        for new_link in new_links:
            if new_link in segment_contigs:
                for contig in segment_contigs[new_link]:
                    if i != contigs_map_rev[int(contig)]:
                        # Add edge to list of edges
                        edge_list.append((i, contigs_map_rev[int(contig)]))

    return edge_list


def get_links_flye(assembly_graph_file):

    my_map = BidirectionalMap()

    node_count = 0

    links = []

    # Get contig connections from .gfa file
    with open(assembly_graph_file) as file:
        line = file.readline()

        while line != "":

            # Count the number of contigs
            if "S" in line:
                strings = line.split("\t")
                my_node = strings[1]
                my_map[node_count] = my_node
                node_count += 1

            # Identify lines with link information
            elif "L" in line:

                link = []
                strings = line.split("\t")

                if strings[1] != strings[3]:
                    start = strings[1]
                    end = strings[3]
                    link.append(start)
                    link.append(end)
                    links.append(link)

            line = file.readline()

    return node_count, links, my_map


def get_graph_edges_flye(links, contig_names_rev):

    edge_list = []

    for link in links:
        # Remove self loops
        if link[0] != link[1]:
            # Add edge to list of edges
            edge_list.append(
                (contig_names_rev[link[0]], contig_names_rev[link[1]]))

    return edge_list


def get_links_megahit(assembly_graph_file):

    node_count = 0

    graph_contigs = {}

    links = []

    my_map = BidirectionalMap()

    # Get links from .gfa file
    with open(assembly_graph_file) as file:

        line = file.readline()

        while line != "":

            # Identify lines with link information
            if line.startswith("L"):
                link = []

                strings = line.split("\t")

                link1 = strings[1]
                link2 = strings[3]

                link.append(link1)
                link.append(link2)
                links.append(link)

            elif line.startswith("S"):
                strings = line.split()

                my_map[node_count] = strings[1]

                graph_contigs[strings[1]] = strings[2]

                node_count += 1

            line = file.readline()

    return node_count, graph_contigs, links, my_map


def get_graph_edges_megahit(links, contig_names_rev):

    edge_list = []

    # Iterate links
    for link in links:
        # Remove self loops
        if link[0] != link[1]:
            # Add edge to list of edges
            edge_list.append(
                (contig_names_rev[link[0]], contig_names_rev[link[1]]))

    return edge_list


def get_isolated(node_count, assembly_graph):

    isolated = []

    # Get isolated contigs which have no neighbours
    for i in range(node_count):

        neighbours = assembly_graph.neighbors(i, mode="ALL")

        if len(neighbours) == 0:
            isolated.append(i)

    return isolated


def get_non_isolated(node_count, assembly_graph, binned_contigs):

    non_isolated = []

    # Get connectes contigs labelled components
    for i in range(node_count):

        if i not in non_isolated and i in binned_contigs:

            component = []
            component.append(i)
            length = len(component)
            neighbours = assembly_graph.neighbors(i, mode="ALL")

            for neighbor in neighbours:
                if neighbor not in component:
                    component.append(neighbor)

            component = list(set(component))

            while length != len(component):

                length = len(component)

                for j in component:

                    neighbours = assembly_graph.neighbors(j, mode="ALL")

                    for neighbor in neighbours:
                        if neighbor not in component:
                            component.append(neighbor)

            labelled = False
            for j in component:
                if j in binned_contigs:
                    labelled = True
                    break

            if labelled:
                for j in component:
                    if j not in non_isolated:
                        non_isolated.append(j)

    return non_isolated
