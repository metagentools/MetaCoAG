#!/usr/bin/env python3

import itertools
import logging
import pickle
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from Bio import SeqIO

__author__ = "Vijini Mallawaarachchi and Yu Lin"
__copyright__ = "Copyright 2020, MetaCoAG Project"
__license__ = "GPL-3.0"
__version__ = "1.2.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Stable Release"


# Create logger
logger = logging.getLogger(f"MetaCoaAG {__version__}")

# Set complements of each nucleotide
complements = {"A": "T", "C": "G", "G": "C", "T": "A"}

# Set bits for each nucleotide
nt_bits = {"A": 0, "C": 1, "G": 2, "T": 3}

VERY_SMALL_VAL = 0.0001


def get_rc(seq):
    rev = reversed(seq)
    return "".join([complements.get(i, i) for i in rev])


def mer2bits(kmer):
    bit_mer = nt_bits.get(kmer[0], 0)
    for c in kmer[1:]:
        bit_mer = (bit_mer << 2) | nt_bits.get(c, 0)
    return bit_mer


def compute_kmer_inds(k):
    kmer_inds = {}
    kmer_count_len = 0

    alphabet = "ACGT"

    all_kmers = ["".join(kmer) for kmer in itertools.product(alphabet, repeat=k)]
    all_kmers.sort()
    ind = 0
    for kmer in all_kmers:
        bit_mer = mer2bits(kmer)
        rc_bit_mer = mer2bits(get_rc(kmer))
        if rc_bit_mer in kmer_inds:
            kmer_inds[bit_mer] = kmer_inds[rc_bit_mer]
        else:
            kmer_inds[bit_mer] = ind
            kmer_count_len += 1
            ind += 1

    return kmer_inds, kmer_count_len


def count_kmers(args):
    contig_num, seq, k, kmer_inds, kmer_count_len = args
    profile = np.zeros(kmer_count_len)
    seq = seq.strip()

    for i in range(0, len(seq) - k + 1):
        bit_mer = mer2bits(seq[i : (i + k)])
        index = kmer_inds[bit_mer]
        profile[index] += 1

    return contig_num, profile / max(1, sum(profile))


def _get_tetramer_cache_path(output_path, contigs_file):
    contigs_name = Path(contigs_file).name
    return Path(output_path) / f"{contigs_name}.normalized_contig_tetramers.pickle"


def _get_contigs_file_metadata(contigs_file):
    contigs_path = Path(contigs_file)
    stat = contigs_path.stat()
    return {
        "path": str(contigs_path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _get_tetramer_cache_signature(contigs_file, min_length, contig_lengths):
    return {
        "contigs_file": _get_contigs_file_metadata(contigs_file),
        "min_length": min_length,
        "eligible_contigs": tuple(
            sorted(
                (contig, length)
                for contig, length in contig_lengths.items()
                if length >= min_length
            )
        ),
    }


def _load_cached_tetramer_profiles(cache_path, signature):
    if not cache_path.is_file():
        return None

    with open(cache_path, "rb") as handle:
        cached = pickle.load(handle)

    if not isinstance(cached, dict) or cached.get("version") != 2:
        return None

    if cached.get("signature") != signature:
        return None

    profiles = cached.get("profiles")
    if not isinstance(profiles, dict):
        return None

    return profiles


def _iter_tetramer_work_items(
    contigs_file,
    contig_names_rev,
    graph_to_contig_map_rev,
    contig_lengths,
    min_length,
    kmer_inds,
    kmer_count_len,
):
    for record in SeqIO.parse(contigs_file, "fasta"):
        if graph_to_contig_map_rev is None:
            if record.id not in contig_names_rev:
                continue
            contig_num = contig_names_rev[record.id]
        else:
            if record.id not in graph_to_contig_map_rev:
                continue
            contig_num = contig_names_rev[graph_to_contig_map_rev[record.id]]

        if contig_lengths[contig_num] >= min_length:
            yield (
                contig_num,
                str(record.seq),
                4,
                kmer_inds,
                kmer_count_len,
            )


def get_tetramer_profiles(
    output_path,
    contigs_file,
    contig_names_rev,
    contig_lengths,
    min_length,
    nthreads,
    graph_to_contig_map_rev=None,
):
    cache_path = _get_tetramer_cache_path(output_path, contigs_file)
    signature = _get_tetramer_cache_signature(contigs_file, min_length, contig_lengths)
    normalized_tetramer_profiles = _load_cached_tetramer_profiles(cache_path, signature)

    if normalized_tetramer_profiles is None:
        normalized_tetramer_profiles = {}
        kmer_inds_4, kmer_count_len_4 = compute_kmer_inds(4)
        work_items = _iter_tetramer_work_items(
            contigs_file=contigs_file,
            contig_names_rev=contig_names_rev,
            graph_to_contig_map_rev=graph_to_contig_map_rev,
            contig_lengths=contig_lengths,
            min_length=min_length,
            kmer_inds=kmer_inds_4,
            kmer_count_len=kmer_count_len_4,
        )

        if nthreads == 1:
            for contig_num, normalized_profile in map(count_kmers, work_items):
                normalized_tetramer_profiles[contig_num] = normalized_profile
        else:
            with Pool(nthreads) as pool:
                for contig_num, normalized_profile in pool.imap_unordered(
                    count_kmers, work_items, chunksize=100
                ):
                    normalized_tetramer_profiles[contig_num] = normalized_profile

        with open(cache_path, "wb") as handle:
            pickle.dump(
                {
                    "version": 2,
                    "signature": signature,
                    "profiles": normalized_tetramer_profiles,
                },
                handle,
                protocol=pickle.HIGHEST_PROTOCOL,
            )

    return normalized_tetramer_profiles


def get_cov_len(contigs_file, contig_names_rev, min_length, abundance_file):
    coverages = defaultdict(lambda: [0])

    contig_lengths = defaultdict(int)

    for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
        contig_num = contig_names_rev[record.id]

        length = len(record.seq)

        contig_lengths[contig_num] = length

    with open(abundance_file, "r") as my_abundance:
        for line in my_abundance:
            strings = line.strip().split("\t")

            contig_num = contig_names_rev[strings[0]]

            if contig_lengths[contig_num] >= min_length:
                for i in range(1, len(strings)):
                    contig_coverage = float(strings[i])

                    if contig_coverage < VERY_SMALL_VAL:
                        contig_coverage = VERY_SMALL_VAL

                    if contig_num not in coverages:
                        coverages[contig_num] = [contig_coverage]
                    else:
                        coverages[contig_num].append(contig_coverage)

    if len(coverages) == 0:
        logger.error(f"Could not find any contigs longer than {min_length}bp.")
        logger.info("Exiting MetaCoAG... Bye...!")
        sys.exit(1)

    sample_vals = list(coverages.keys())
    n_samples = len(coverages[sample_vals[0]])

    return coverages, contig_lengths, n_samples


def get_cov_len_megahit(
    contigs_file, contig_names_rev, graph_to_contig_map_rev, min_length, abundance_file
):
    coverages = {}

    contig_lengths = {}

    for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
        if record.id not in graph_to_contig_map_rev:
            continue
        contig_num = contig_names_rev[graph_to_contig_map_rev[record.id]]
        length = len(record.seq)
        contig_lengths[contig_num] = length

    with open(abundance_file, "r") as my_abundance:
        for line in my_abundance:
            strings = line.strip().split("\t")

            if strings[0] not in graph_to_contig_map_rev:
                continue
            contig_num = contig_names_rev[graph_to_contig_map_rev[strings[0]]]

            if contig_lengths[contig_num] >= min_length:
                for i in range(1, len(strings)):
                    contig_coverage = float(strings[i])

                    if contig_coverage < VERY_SMALL_VAL:
                        contig_coverage = VERY_SMALL_VAL

                    if contig_num not in coverages:
                        coverages[contig_num] = [contig_coverage]
                    else:
                        coverages[contig_num].append(contig_coverage)

    if len(coverages) == 0:
        logger.error(f"Could not find any contigs longer than {min_length}bp.")
        logger.info("Exiting MetaCoAG... Bye...!")
        sys.exit(1)

    n_samples = len(coverages[list(coverages.keys())[0]])

    return coverages, contig_lengths, n_samples


def get_bin_profiles(bins, coverages, normalized_tetramer_profiles):
    bin_tetramer_profile = {}
    bin_coverage_profile = {}

    for b in bins:
        coverage_b = []
        tetramer_b = []

        for contig in bins[b]:
            coverage_b.append(coverages[contig])
            tetramer_b.append(normalized_tetramer_profiles[contig])

        coverage_b = np.array(coverage_b)
        tetramer_b = np.array(tetramer_b)

        bin_coverage_profile[b] = np.mean(coverage_b, axis=0)
        bin_tetramer_profile[b] = np.mean(tetramer_b, axis=0)

    return bin_tetramer_profile, bin_coverage_profile
