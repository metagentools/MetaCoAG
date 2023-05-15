import os
import subprocess
from pathlib import Path

import pytest

__author__ = "Vijini Mallawaarachchi and Yu Lin"
__copyright__ = "Copyright 2020, MetaCoAG Project"
__license__ = "GPL-3.0"
__version__ = "1.1.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Stable Release"


TEST_ROOTDIR = Path(__file__).parent
EXEC_ROOTDIR = Path(__file__).parent.parent


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def get_files_and_seq_counts(output_path):
    output_files = os.listdir(output_path)
    seq_counts = []
    for file in output_files:
        seq_count = 0
        with open(f"{output_path}/{file}", "r") as myfile:
            for line in myfile:
                if line.strip().startswith(">"):
                    seq_count += 1
        seq_counts.append(seq_count)

    seq_counts.sort()

    return len(output_files), seq_counts


@pytest.fixture(scope="session")
def test_metacoag_spades_command(tmp_path_factory):
    """test metacoag on metaspades assembly"""

    tmp_dir = tmp_path_factory.mktemp("tmp")

    dir_name = TEST_ROOTDIR / "data" / "5G_metaspades"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    abundance = dir_name / "coverm_mean_coverage.tsv"
    cmd = f"metacoag --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --abundance {abundance} --output {tmp_dir}"
    exec_command(cmd)

    n_bins, seq_counts = get_files_and_seq_counts(tmp_dir / "bins")

    return n_bins, seq_counts


def test_n_bins_metacoag_spades(test_metacoag_spades_command):
    n_bins, seq_counts = test_metacoag_spades_command

    # Assert number of bins
    assert n_bins == 5

    # Assert bin sizes
    assert seq_counts == [10, 23, 48, 69, 78]


@pytest.fixture(scope="session")
def test_metacoag_megahit_command(tmp_path_factory):
    """test metacoag on megahit assembly"""

    tmp_dir = tmp_path_factory.mktemp("tmp")

    dir_name = TEST_ROOTDIR / "data" / "5G_MEGAHIT"
    graph = dir_name / "final.gfa"
    contigs = dir_name / "final.contigs.fa"
    abundance = dir_name / "abundance.tsv"
    cmd = f"metacoag --assembler megahit --graph {graph} --contigs {contigs} --abundance {abundance} --output {tmp_dir}"
    exec_command(cmd)

    n_bins, seq_counts = get_files_and_seq_counts(tmp_dir / "bins")

    return n_bins, seq_counts


def test_n_bins_metacoag_megahit(test_metacoag_megahit_command):
    n_bins, seq_counts = test_metacoag_megahit_command

    # Assert number of bins
    assert n_bins == 5

    # Assert bin sizes
    assert seq_counts == [36, 40, 46, 84, 127]
