#!/usr/bin/env python3

import logging
import pathlib
import subprocess

__author__ = "Vijini Mallawaarachchi and Yu Lin"
__copyright__ = "Copyright 2020, MetaCoAG Project"
__license__ = "GPL-3.0"
__version__ = "1.2.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Stable Release"


# create logger
logger = logging.getLogger(f"MetaCoaAG {__version__}")


# Modified from SolidBin
def scan_for_marker_genes(contigs_file, nthreads, markerURL, no_cut_tc, hard=0):
    if markerURL == "auxiliary/marker.hmm":
        markerURL = pathlib.Path(__file__).parent / "auxiliary" / "marker.hmm"

    marker_path = pathlib.Path(markerURL)
    contigs_path = pathlib.Path(contigs_file)
    frag_prefix = pathlib.Path(f"{contigs_file}.frag")
    frag_result = pathlib.Path(f"{frag_prefix}.faa")
    hmm_result = pathlib.Path(f"{contigs_file}.hmmout")

    logger.info("Using marker file: %s", marker_path)

    if not frag_result.exists():
        frag_command = [
            "run_FragGeneScan.pl",
            f"-genome={contigs_path}",
            f"-out={frag_prefix}",
            "-complete=0",
            "-train=complete",
            f"-thread={nthreads}",
        ]
        logger.debug("exec cmd: %s", subprocess.list2cmdline(frag_command))
        with open(f"{frag_prefix}.out", "w") as stdout, open(
            f"{frag_prefix}.err", "w"
        ) as stderr:
            subprocess.run(frag_command, stdout=stdout, stderr=stderr, check=True)

    if not frag_result.exists():
        raise FileNotFoundError(
            f"FragGeneScan completed without creating expected output: {frag_result}"
        )

    if hmm_result.exists():
        return

    hmm_command = ["hmmsearch", "--domtblout", str(hmm_result)]
    if not no_cut_tc:
        hmm_command.append("--cut_tc")
    hmm_command.extend(
        ["--cpu", str(nthreads), str(marker_path), str(frag_result)]
    )
    logger.debug("exec cmd: %s", subprocess.list2cmdline(hmm_command))
    with open(f"{hmm_result}.out", "w") as stdout, open(
        f"{hmm_result}.err", "w"
    ) as stderr:
        subprocess.run(hmm_command, stdout=stdout, stderr=stderr, check=True)

    if not hmm_result.exists():
        raise FileNotFoundError(
            f"hmmsearch completed without creating expected output: {hmm_result}"
        )


# Get contigs containing marker genes
def get_contigs_with_marker_genes(
    contigs_file, contig_names_rev, mg_length_threshold, contig_lengths, min_length
):
    marker_contigs = {}
    marker_contig_counts = {}
    contig_markers = {}

    with open(f"{contigs_file}.hmmout", "r") as myfile:
        for line in myfile:
            if not line.startswith("#"):
                strings = line.strip().split()

                contig = strings[0]

                # Marker gene name
                marker_gene = strings[3]

                # Marker gene length
                marker_gene_length = int(strings[5])

                # Mapped marker gene length
                mapped_marker_length = int(strings[16]) - int(strings[15])

                name_strings = contig.split("_")
                name_strings = name_strings[: len(name_strings) - 3]

                # Contig name
                contig_name = "_".join(name_strings)

                contig_num = contig_names_rev[contig_name]
                contig_length = contig_lengths[contig_num]

                if (
                    contig_length >= min_length
                    and mapped_marker_length > marker_gene_length * mg_length_threshold
                ):
                    marker_repeated_in_contig = False

                    contig_markers.setdefault(contig_num, set()).add(marker_gene)

                    # Get contigs containing each marker gene
                    if marker_gene not in marker_contigs:
                        marker_contigs[marker_gene] = [contig_num]
                    else:
                        if contig_num not in marker_contigs[marker_gene]:
                            marker_contigs[marker_gene].append(contig_num)
                        else:
                            marker_repeated_in_contig = True

                    # Get contig counts for each marker
                    if marker_gene not in marker_contig_counts:
                        marker_contig_counts[marker_gene] = 1
                    else:
                        if not marker_repeated_in_contig:
                            marker_contig_counts[marker_gene] += 1

    return marker_contigs, marker_contig_counts, contig_markers


# Get contigs containing marker genes
def get_contigs_with_marker_genes_megahit(
    contigs_file,
    contig_names_rev,
    graph_to_contig_map_rev,
    mg_length_threshold,
    contig_lengths,
    min_length,
):
    marker_contigs = {}
    marker_contig_counts = {}
    contig_markers = {}

    with open(f"{contigs_file}.hmmout", "r") as myfile:
        for line in myfile:
            if not line.startswith("#"):
                strings = line.strip().split()

                contig = strings[0]

                # Marker gene name
                marker_gene = strings[3]

                # Marker gene length
                marker_gene_length = int(strings[5])

                # Mapped marker gene length
                mapped_marker_length = int(strings[16]) - int(strings[15])

                name_strings = contig.split("_")
                name_strings = name_strings[: len(name_strings) - 3]

                # Contig name
                contig_name = "_".join(name_strings)

                if contig_name not in graph_to_contig_map_rev:
                    continue
                contig_num = contig_names_rev[graph_to_contig_map_rev[contig_name]]
                contig_length = contig_lengths[contig_num]

                if (
                    contig_length >= min_length
                    and mapped_marker_length > marker_gene_length * mg_length_threshold
                ):
                    marker_repeated_in_contig = False

                    contig_markers.setdefault(contig_num, set()).add(marker_gene)

                    # Get contigs containing each marker gene
                    if marker_gene not in marker_contigs:
                        marker_contigs[marker_gene] = [contig_num]
                    else:
                        if contig_num not in marker_contigs[marker_gene]:
                            marker_contigs[marker_gene].append(contig_num)
                        else:
                            marker_repeated_in_contig = True

                    # Get contig counts for each marker
                    if marker_gene not in marker_contig_counts:
                        marker_contig_counts[marker_gene] = 1
                    else:
                        if not marker_repeated_in_contig:
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
