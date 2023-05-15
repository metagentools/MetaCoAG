#!/usr/bin/env python3

import multiprocessing as mp
import re
from collections import defaultdict

from Bio import SeqIO

from metacoag_utils.bidirectionalmap import BidirectionalMap


def get_segment_paths_spades(contig_paths):
    paths = {}
    segment_contigs = {}
    node_count = 0

    my_map = BidirectionalMap()

    contig_names = BidirectionalMap()

    current_contig_num = ""

    with open(contig_paths) as file:
        name = file.readline().strip()
        path = file.readline().strip()

        while name != "" and path != "":
            while ";" in path:
                path = path[:-2] + "," + file.readline()

            start = "NODE_"
            end = "_length_"
            contig_num = str(int(re.search("%s(.*)%s" % (start, end), name).group(1)))

            segments = path.rstrip().split(",")

            if current_contig_num != contig_num:
                my_map[node_count] = int(contig_num)
                contig_names[node_count] = name.strip()
                current_contig_num = contig_num
                node_count += 1

            if contig_num not in paths:
                paths[contig_num] = segments

            for segment in segments:
                if segment not in segment_contigs:
                    segment_contigs[segment] = set([contig_num])
                else:
                    segment_contigs[segment].add(contig_num)

            name = file.readline().strip()
            path = file.readline().strip()

    return paths, segment_contigs, node_count, my_map, contig_names


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
                strings = line.split("\t")
                f1, f2 = strings[1] + strings[2], strings[3] + strings[4]
                links_map[f1].add(f2)
                links_map[f2].add(f1)
                links.append(strings[1] + strings[2] + " " + strings[3] + strings[4])
            line = file.readline()

    # Create list of edges
    edge_list = []

    for i in range(len(paths)):
        segments = paths[str(contigs_map[i])]

        new_links = []

        for segment in segments:
            my_segment = segment

            my_segment_rev = ""

            if my_segment.endswith("+"):
                my_segment_rev = my_segment[:-1] + "-"
            else:
                my_segment_rev = my_segment[:-1] + "+"

            if segment in links_map:
                new_links.extend(list(links_map[segment]))

            if my_segment_rev in links_map:
                new_links.extend(list(links_map[my_segment_rev]))

        if my_segment in segment_contigs:
            for contig in segment_contigs[my_segment]:
                if i != contigs_map_rev[int(contig)]:
                    # Add edge to list of edges
                    edge_list.append((i, contigs_map_rev[int(contig)]))

        if my_segment_rev in segment_contigs:
            for contig in segment_contigs[my_segment_rev]:
                if i != contigs_map_rev[int(contig)]:
                    # Add edge to list of edges
                    edge_list.append((i, contigs_map_rev[int(contig)]))

        for new_link in new_links:
            if new_link in segment_contigs:
                for contig in segment_contigs[new_link]:
                    if i != contigs_map_rev[int(contig)]:
                        # Add edge to list of edges
                        edge_list.append((i, contigs_map_rev[int(contig)]))

    return edge_list


def get_flye_contig_map(contigs_file):
    contig_names = BidirectionalMap()

    contig_num = 0

    for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
        contig_names[contig_num] = record.id
        contig_num += 1

    return contig_names


def get_links_flye(contig_paths, contig_names_rev):
    paths = {}
    segment_contigs = {}

    my_map = BidirectionalMap()

    with open(contig_paths) as file:
        for line in file.readlines():
            if not line.startswith("#"):
                strings = line.strip().split()

                contig_name = strings[0]

                path = strings[-1]
                path = path.replace("*", "")

                if path.startswith(","):
                    path = path[1:]

                if path.endswith(","):
                    path = path[:-1]

                segments = path.rstrip().split(",")

                contig_num = contig_names_rev[contig_name]

                if contig_num not in paths:
                    paths[contig_num] = segments

                for segment in segments:
                    if segment not in segment_contigs:
                        segment_contigs[segment] = set([contig_num])
                    else:
                        segment_contigs[segment].add(contig_num)

    return paths, segment_contigs, len(contig_names_rev), my_map


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


def get_links_megahit_custom(assembly_graph_file):
    my_map = BidirectionalMap()

    node_count = 0

    nodes = []

    links = []

    # Get contig connections from .gfa file
    with open(assembly_graph_file) as file:
        for line in file.readlines():
            line = line.strip()

            # Count the number of contigs
            if line.startswith("S"):
                strings = line.split("\t")
                my_node = strings[1][:-2]
                my_map[node_count] = my_node
                nodes.append(my_node)
                node_count += 1

            # Identify lines with link information
            elif line.startswith("L"):
                link = []
                strings = line.split("\t")

                if strings[1] != strings[3]:
                    start = strings[1]
                    end = strings[3]
                    link.append(start)
                    link.append(end)
                    links.append(link)

    return node_count, links, my_map


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


def get_connected_components(i, assembly_graph, binned_contigs):
    non_isolated = []
    if i not in non_isolated and i in binned_contigs:
        component = []
        component.append(i)
        length = len(component)
        neighbours = assembly_graph.neighbors(i, mode="ALL")
        for neighbour in neighbours:
            if neighbour not in component:
                component.append(neighbour)
        component = list(set(component))
        while length != len(component):
            length = len(component)
            for j in component:
                neighbours = assembly_graph.neighbors(j, mode="ALL")
                for neighbour in neighbours:
                    if neighbour not in component:
                        component.append(neighbour)
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


def get_non_isolated(node_count, assembly_graph, binned_contigs, nthreads):
    with mp.Pool(processes=nthreads) as pool:
        non_isolated = pool.starmap(
            get_connected_components,
            [(i, assembly_graph, binned_contigs) for i in range(node_count)],
        )
    return non_isolated
