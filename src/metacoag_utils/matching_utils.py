import math

from scipy.spatial import distance
from scipy.stats import poisson

mu_intra, sigma_intra = 0, 0.01037897 / 2
mu_inter, sigma_inter = 0.0676654, 0.03419337


def normpdf(x, mean, sd):
    var = float(sd)**2
    denom = sd*(2*math.pi)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom


def get_tetramer_distance(seq1, seq2):
    return distance.euclidean(seq1, seq2)


def get_comp_probability(tetramer_dist):
    return normpdf(tetramer_dist, mu_intra, sigma_intra)/(normpdf(tetramer_dist, mu_intra, sigma_intra)+normpdf(tetramer_dist, mu_inter, sigma_inter))


def get_cov_probability(cov1, cov2):
    return poisson.pmf(cov1, cov2)