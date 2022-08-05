import subprocess
from pathlib import Path

import pytest

__author__ = "Vijini Mallawaarachchi and Yu Lin"
__copyright__ = "Copyright 2020, MetaCoAG Project"
__license__ = "GPL-3.0"
__version__ = "1.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Development"


TEST_ROOTDIR = Path(__file__).parent
EXEC_ROOTDIR = Path(__file__).parent.parent


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    """set the working directory for all tests"""
    monkeypatch.chdir(tmp_dir)


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


def test_metacoag_spades_command(tmp_dir):
    """test metacoag on spades assembly"""
    dir_name = TEST_ROOTDIR / "data" / "5G_metaspades"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    abundance = dir_name / "coverm_mean_coverage.tsv"
    cmd = f"{EXEC_ROOTDIR}/metacoag --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --abundance {abundance} --output {tmp_dir}"
    exec_command(cmd)


def test_metacoag_megahit_command(tmp_dir):
    """test metacoag on megahit assembly"""
    dir_name = TEST_ROOTDIR / "data" / "5G_MEGAHIT"
    graph = dir_name / "final.gfa"
    contigs = dir_name / "final.contigs.fa"
    abundance = dir_name / "abundance.tsv"
    cmd = f"{EXEC_ROOTDIR}/metacoag --assembler megahit --graph {graph} --contigs {contigs} --abundance {abundance} --output {tmp_dir}"
    exec_command(cmd)
