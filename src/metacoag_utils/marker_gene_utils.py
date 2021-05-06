#!/usr/bin/env python3

import os
import pathlib
import re
import numpy as np
import logging

# create logger
logger = logging.getLogger('MetaCoaAG 0.1')


# Modified from SolidBin
def scan_for_marker_genes(contig_file, nthreads, hard=0):

    software_path = pathlib.Path(__file__).parent.parent.absolute()

    fragScanURL = os.path.join(software_path.parent, 'auxiliary',
                               'FragGeneScan1.31', 'run_FragGeneScan.pl')
    hmmExeURL = os.path.join(software_path.parent, 'auxiliary', 'hmmer-3.3',
                             'src', 'hmmsearch')
    markerURL = os.path.join(software_path.parent, 'auxiliary', 'marker.hmm')

    logger.debug(markerURL)

    fragResultURL = contig_file+".frag.faa"
    hmmResultURL = contig_file+".hmmout"
    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL+" -genome="+contig_file+" -out="+contig_file + \
            ".frag -complete=0 -train=complete -thread="+str(nthreads)+" 1>" + \
            contig_file+".frag.out 2>"+contig_file+".frag.err"
        logger.debug("exec cmd: "+fragCmd)
        os.system(fragCmd)

    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL+" --domtblout "+hmmResultURL+" --cut_tc --cpu "+str(nthreads)+" " + \
                markerURL+" "+fragResultURL+" 1>"+hmmResultURL+".out 2>"+hmmResultURL+".err"
            logger.debug("exec cmd: "+hmmCmd)
            os.system(hmmCmd)

        else:
            logger.debug("HMMER search failed! Path: " +
                         hmmResultURL + " does not exist.")
    else:
        logger.debug("FragGeneScan failed! Path: " +
                     fragResultURL + " does not exist.")


def get_contigs_with_marker_genes(contigs_file, contig_names_rev, mg_length_threshold, contig_lengths, min_length):

    marker_contigs = {}
    marker_contig_counts = {}
    contig_markers = {}

    with open(contigs_file+".hmmout", "r") as myfile:
        for line in myfile.readlines():
            if not line.startswith("#"):
                strings = line.strip().split()

                contig = strings[0]
                marker_gene = strings[3]
                marker_gene_length = int(strings[5])

                mapped_marker_length = int(strings[16]) - int(strings[15])

                name_strings = contig.split("_")
                name_strings = name_strings[:len(name_strings)-3]

                contig_name = "_".join(name_strings)

                contig_num = contig_names_rev[contig_name]
                contig_length = contig_lengths[contig_num]

                if contig_length >= min_length and mapped_marker_length > marker_gene_length*mg_length_threshold:

                    if contig_num not in contig_markers:
                        contig_markers[contig_num] = [marker_gene]
                    else:
                        contig_markers[contig_num].append(marker_gene)

                    if marker_gene not in marker_contigs:
                        marker_contigs[marker_gene] = [contig_num]
                    else:
                        marker_contigs[marker_gene].append(contig_num)

                    if marker_gene not in marker_contig_counts:
                        marker_contig_counts[marker_gene] = 1
                    else:
                        marker_contig_counts[marker_gene] += 1

    return marker_contigs, marker_contig_counts, contig_markers


def count_contigs_with_marker_genes(marker_contig_counts):

    marker_frequencies = {}

    for marker in marker_contig_counts:

        if marker_contig_counts[marker] not in marker_frequencies:
            marker_frequencies[marker_contig_counts[marker]] = 1
        else:
            marker_frequencies[marker_contig_counts[marker]] += 1

    return marker_frequencies


def get_seed_marker_gene_counts(marker_contig_counts, seed_mg_threshold):

    my_gene_counts = []

    vals = list(marker_contig_counts.values())

    data = np.array(vals)

    values, counts = np.unique(data, return_counts=True)

    max_count_gene = values[np.where(counts == np.amax(counts))]

    my_gene_counts.append(max_count_gene[0])

    max_count_gene_index = np.where(values == max_count_gene[0])[0][0]

    min_count_index = max_count_gene_index
    max_count_index = max_count_gene_index

    for i in range(max_count_gene_index+1, len(values)):
        if counts[max_count_gene_index]*seed_mg_threshold <= counts[i]:
            max_count_index = i

    for i in range(max_count_gene_index, 0, -1):
        if counts[max_count_gene_index]*seed_mg_threshold <= counts[i-1]:
            min_count_index = i-1

    for i in range(min_count_index, max_count_index+1):
        if values[i] not in my_gene_counts:
            my_gene_counts.append(values[i])

    return my_gene_counts
