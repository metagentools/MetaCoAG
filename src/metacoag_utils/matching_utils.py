#!/usr/bin/env python3

import math

from scipy.spatial import distance

MU_INTRA, SIGMA_INTRA = 0, 0.01037897 / 2
MU_INTER, SIGMA_INTER = 0.0676654, 0.03419337
VERY_SMALL_DOUBLE = 1e-100


def normpdf(x, mean, sd):
    var = float(sd)**2
    denom = sd*(2*math.pi)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom


def get_tetramer_distance(seq1, seq2):
    return distance.euclidean(seq1, seq2)


def get_comp_probability(tetramer_dist):
    gaus_intra = normpdf(tetramer_dist, MU_INTRA, SIGMA_INTRA)
    gaus_inter = normpdf(tetramer_dist, MU_INTER, SIGMA_INTER)
    return gaus_intra/(gaus_intra+gaus_inter)


def get_cov_probability(cov1, cov2):

    poisson_prod = 1

    for i in range(len(cov1)):

        if cov1[i] != 0.0 and cov2[i] != 0.0:

            # Adapted from http://www.masaers.com/2013/10/08/Implementing-Poisson-pmf.html
            poisson_pmf = math.exp(
                (cov1[i] * math.log(cov2[i])) - math.lgamma(cov1[i] + 1.0) - cov2[i])

            if poisson_pmf < VERY_SMALL_DOUBLE:
                poisson_pmf = VERY_SMALL_DOUBLE

            poisson_prod = poisson_prod * poisson_pmf

    if poisson_prod == 1:
        poisson_prod = VERY_SMALL_DOUBLE

    return poisson_prod
