#!/usr/bin/env python3

import numpy as np
import itertools

complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
nt_bits = {'A':0,'C':1,'G':2,'T':3}

def get_rc(seq):
    rev = reversed(seq)
    return "".join([complements.get(i,i) for i in rev])


def mer2bits(kmer):
    bit_mer=nt_bits[kmer[0]]
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
    profile = profile/max(1, sum(profile))
    
    return profile