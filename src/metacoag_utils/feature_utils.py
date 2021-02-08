#!/usr/bin/env python3

import numpy as np
import re
import itertools
import os

from Bio import SeqIO
from multiprocessing import Pool

complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
nt_bits = {'A':0,'C':1,'G':2,'T':3}

def get_rc(seq):
    rev = reversed(seq)
    return "".join([complements.get(i,i) for i in rev])


def mer2bits(kmer):
    bit_mer = nt_bits[kmer[0]]
    for c in kmer[1:]:
        bit_mer = (bit_mer << 2) | nt_bits[c]
    return bit_mer

def compute_kmer_inds(k):
    kmer_inds = {}
    kmer_count_len = 0

    alphabet = 'ACGT'
    
    all_kmers = [''.join(kmer) for kmer in itertools.product(alphabet,repeat=k)]
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
    seq, k, kmer_inds, kmer_count_len = args
    profile = np.zeros(kmer_count_len)
    arrs = []
    seq = list(seq.strip())
    
    for i in range(0, len(seq) - k + 1):
        bit_mer = mer2bits(seq[i:i+k])
        index = kmer_inds[bit_mer]
        profile[index] += 1
    # profile = profile/max(1, sum(profile))
    
    return profile, profile/max(1, sum(profile))


def get_tetramer_profiles(output_path, seqs, nthreads):

    tetramer_profiles = {}
    normalized_tetramer_profiles = {}

    if os.path.isfile(output_path+"contig_tetramers.txt") and 
            os.path.isfile(output_path+"normalized_contig_tetramers.txt"):
        i = 0
        with open(output_path+"contig_tetramers.txt") as tetramers_file:
            for line in tetramers_file.readlines():
                f_list = [float(i) for i in line.split(" ") if i.strip()]
                tetramer_profiles[i] = f_list
                i += 1

        i = 0
        with open(output_path+"normalized_contig_tetramers.txt") as tetramers_file:
            for line in tetramers_file.readlines():
                f_list = [float(i) for i in line.split(" ") if i.strip()]
                normalized_tetramer_profiles[i] = f_list
                i += 1

    else:

        kmer_inds_4, kmer_count_len_4 = compute_kmer_inds(4)

        pool = Pool(nthreads)
        record_tetramers = pool.map(count_kmers, 
                                    [(seq, 4, kmer_inds_4, kmer_count_len_4) for seq in seqs])
        pool.close()

        normalized = [x[1] for x in record_tetramers]
        unnormalized = [x[0] for x in record_tetramers]
        
        i = 0

        for l in range(len(unnormalized)):
            tetramer_profiles[i] = unnormalized[l]
            i += 1
        
        with open(output_path+"contig_tetramers.txt", "w+") as myfile:
            for l in range(len(unnormalized)):
                for j in range(len(unnormalized[l])):
                    myfile.write(str(unnormalized[l][j])+" ")
                myfile.write("\n")

        i = 0

        for l in range(len(normalized)):
            normalized_tetramer_profiles[i] = normalized[l]
            i += 1
        
        with open(output_path+"normalized_contig_tetramers.txt", "w+") as myfile:
            for l in range(len(normalized)):
                for j in range(len(normalized[l])):
                    myfile.write(str(normalized[l][j])+" ")
                myfile.write("\n")

    return tetramer_profiles, normalized_tetramer_profiles


def get_cov_len_spades(contigs_file, contigs_map_rev):

    coverages = {}

    contig_lengths = {}

    i = 0

    seqs = []

    for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
        
        start = 'NODE_'
        end = '_length_'
        contig_num = contigs_map_rev[int(re.search('%s(.*)%s' % (start, end), record.id).group(1))]
        
        start = '_cov_'
        end = ''
        coverage = int(float(re.search('%s(.*)%s' % (start, end), record.id).group(1)))
        
        start = '_length_'
        end = '_cov'
        length = int(re.search('%s(.*)%s' % (start, end), record.id).group(1))
        
        coverages[contig_num] = coverage
        contig_lengths[contig_num] = length
        
        seqs.append(str(record.seq))
        
        i += 1
    
    return seqs, coverages, contig_lengths