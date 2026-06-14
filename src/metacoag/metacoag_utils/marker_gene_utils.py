#!/usr/bin/env python3

import logging
import subprocess

from pathlib import Path


__author__ = "Vijini Mallawaarachchi and Yu Lin"
__copyright__ = "Copyright 2020, MetaCoAG Project"
__license__ = "GPL-3.0"
__version__ = "1.2.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Stable Release"


logger = logging.getLogger(f"MetaCoAG {__version__}")


# Modified from SolidBin
def scan_for_marker_genes(
    contigs_file,
    nthreads,
    marker_url=None,
    no_cut_tc=False,
    hard=0,
    **legacy_options,
):
    """Run FragGeneScan and HMMER for the supplied contigs.

    ``markerURL`` remains accepted as a compatibility keyword.
    """
    if marker_url is None:
        marker_url = legacy_options.pop("markerURL", None)
    if legacy_options:
        unexpected = next(iter(legacy_options))
        raise TypeError(f"Unexpected keyword argument: {unexpected}")
    if marker_url is None:
        raise TypeError("marker_url is required")

    del hard

    if marker_url == "auxiliary/marker.hmm":
        marker_url = Path(__file__).parent / "auxiliary" / "marker.hmm"

    marker_path = Path(marker_url)
    contigs_path = Path(contigs_file)
    frag_prefix = Path(f"{contigs_file}.frag")
    frag_result = Path(f"{frag_prefix}.faa")
    hmm_result = Path(f"{contigs_file}.hmmout")

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
        with (
            open(f"{frag_prefix}.out", "w") as stdout,
            open(f"{frag_prefix}.err", "w") as stderr,
        ):
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
    hmm_command.extend(["--cpu", str(nthreads), str(marker_path), str(frag_result)])
    logger.debug("exec cmd: %s", subprocess.list2cmdline(hmm_command))
    with (
        open(f"{hmm_result}.out", "w") as stdout,
        open(f"{hmm_result}.err", "w") as stderr,
    ):
        subprocess.run(hmm_command, stdout=stdout, stderr=stderr, check=True)

    if not hmm_result.exists():
        raise FileNotFoundError(
            f"hmmsearch completed without creating expected output: {hmm_result}"
        )


def _iter_marker_hits(contigs_file):
    with open(f"{contigs_file}.hmmout", "r") as hmm_output:
        for line in hmm_output:
            if line.startswith("#"):
                continue

            fields = line.strip().split()
            protein_name = fields[0]
            marker_gene = fields[3]
            marker_gene_length = int(fields[5])
            mapped_marker_length = int(fields[16]) - int(fields[15])
            contig_name = "_".join(protein_name.split("_")[:-3])

            yield (
                contig_name,
                marker_gene,
                marker_gene_length,
                mapped_marker_length,
            )


def _add_marker_hit(
    marker_gene,
    contig_id,
    marker_contigs,
    marker_contig_counts,
    contig_markers,
):
    contig_markers.setdefault(contig_id, set()).add(marker_gene)
    marker_contig_list = marker_contigs.setdefault(marker_gene, [])

    if contig_id in marker_contig_list:
        return

    marker_contig_list.append(contig_id)
    marker_contig_counts[marker_gene] = marker_contig_counts.get(marker_gene, 0) + 1


def get_contigs_with_marker_genes(
    contigs_file, contig_names_rev, mg_length_threshold, contig_lengths, min_length
):
    marker_contigs = {}
    marker_contig_counts = {}
    contig_markers = {}

    for (
        contig_name,
        marker_gene,
        marker_gene_length,
        mapped_marker_length,
    ) in _iter_marker_hits(contigs_file):
        contig_id = contig_names_rev[contig_name]
        passes_length = contig_lengths[contig_id] >= min_length
        passes_marker_threshold = (
            mapped_marker_length > marker_gene_length * mg_length_threshold
        )
        if passes_length and passes_marker_threshold:
            _add_marker_hit(
                marker_gene,
                contig_id,
                marker_contigs,
                marker_contig_counts,
                contig_markers,
            )

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

    for (
        contig_name,
        marker_gene,
        marker_gene_length,
        mapped_marker_length,
    ) in _iter_marker_hits(contigs_file):
        if contig_name not in graph_to_contig_map_rev:
            continue

        graph_contig_name = graph_to_contig_map_rev[contig_name]
        contig_id = contig_names_rev[graph_contig_name]
        passes_length = contig_lengths[contig_id] >= min_length
        passes_marker_threshold = (
            mapped_marker_length > marker_gene_length * mg_length_threshold
        )
        if passes_length and passes_marker_threshold:
            _add_marker_hit(
                marker_gene,
                contig_id,
                marker_contigs,
                marker_contig_counts,
                contig_markers,
            )

    return marker_contigs, marker_contig_counts, contig_markers


def count_contigs_with_marker_genes(marker_contig_counts):
    marker_frequencies = {}

    for contig_count in marker_contig_counts.values():
        marker_frequencies[contig_count] = marker_frequencies.get(contig_count, 0) + 1

    return marker_frequencies
