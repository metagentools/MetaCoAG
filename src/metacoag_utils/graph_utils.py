#!/usr/bin/env python3

import re

from collections import defaultdict
from metacoag_utils.bidirectionalmap import BidirectionalMap

def get_segment_paths_spades(contig_paths):

    paths = {}
    segment_contigs = {}
    node_count = 0

    my_map = BidirectionalMap()

    contig_names = {}

    current_contig_num = ""

    with open(contig_paths) as file:
        name = file.readline()
        path = file.readline()

        while name != "" and path != "":

            while ";" in path:
                path = path[:-2]+","+file.readline()

            start = 'NODE_'
            end = '_length_'
            contig_num = str(int(re.search('%s(.*)%s' % (start, end), name).group(1)))

            segments = path.rstrip().split(",")

            if current_contig_num != contig_num:
                my_map[node_count] = int(contig_num)
                contig_names[node_count] = name.strip()
                current_contig_num = contig_num
                node_count += 1

            if contig_num not in paths:
                paths[contig_num] = [segments[0], segments[-1]]

            for segment in segments:
                if segment not in segment_contigs:
                    segment_contigs[segment] = set([contig_num])
                else:
                    segment_contigs[segment].add(contig_num)

            name = file.readline()
            path = file.readline()

    return paths, segment_contigs, node_count, my_map, contig_names



def get_graph_edges_spades(
        assembly_graph_file, node_count, contigs_map, 
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
                f1, f2 = strings[1]+strings[2], strings[3]+strings[4]
                links_map[f1].add(f2)
                links_map[f2].add(f1)
                links.append(strings[1]+strings[2]+" "+strings[3]+strings[4])
            line = file.readline()

    # Create list of edges
    edge_list = []

    for i in range(len(paths)):
        segments = paths[str(contigs_map[i])]

        start = segments[0]
        start_rev = ""

        if start.endswith("+"):
            start_rev = start[:-1]+"-"
        else:
            start_rev = start[:-1]+"+"

        end = segments[1]
        end_rev = ""

        if end.endswith("+"):
            end_rev = end[:-1]+"-"
        else:
            end_rev = end[:-1]+"+"

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
                        edge_list.append((i,contigs_map_rev[int(contig)]))

    return edge_list



def get_isolated(node_count, assembly_graph):
    
    isolated=[]

    for i in range(node_count):
        
        neighbours = assembly_graph.neighbors(i, mode="ALL")
        
        if len(neighbours)==0:
            isolated.append(i)

    return isolated


def get_non_isolated(node_count, assembly_graph, binned_contigs):

    non_isolated = []

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

            while length!= len(component):

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