import numpy as np

from metacoag.metacoag_utils import matching_utils


def test_match_contigs_handles_single_seed_iteration():
    result = matching_utils.match_contigs(
        smg_iteration={0: [0]},
        bins={0: [0]},
        n_bins=0,
        bin_of_contig={0: 0},
        binned_contigs_with_markers=[0],
        bin_markers={0: {"marker_1"}},
        contig_markers={0: {"marker_1"}},
        contig_lengths=np.array([2000]),
        contig_names={0: "contig_0"},
        normalized_tetramer_profiles=np.zeros((1, 136)),
        coverages=np.ones((1, 1)),
        assembly_graph=None,
        w_intra=2.0,
        w_inter=4.0,
        d_limit=20,
    )

    bins, bin_of_contig, n_bins, bin_markers, binned_contigs = result

    assert bins == {0: [0]}
    assert bin_of_contig == {0: 0}
    assert n_bins == 1
    assert bin_markers == {0: {"marker_1"}}
    assert binned_contigs == [0]
