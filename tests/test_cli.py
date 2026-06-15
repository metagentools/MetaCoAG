import pathlib
import pickle
from types import SimpleNamespace

import pytest

from click.testing import CliRunner

from metacoag import metacoag_pipeline, metacoag_runner
from metacoag.cli import main
from metacoag.metacoag_utils.bidirectionalmap import BidirectionalMap


__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]


DATADIR = pathlib.Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    """set the working directory for all tests"""
    monkeypatch.chdir(tmp_dir)


@pytest.fixture(scope="session")
def runner():
    """exportrc works correctly."""
    return CliRunner()


def test_metacoag_spades_run(runner, tmp_dir):
    outpath = tmp_dir
    dir_name = DATADIR / "5G_metaspades"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    abundance = dir_name / "coverm_mean_coverage.tsv"
    args = f"--assembler spades --graph {graph} --contigs {contigs} --paths {paths} --abundance {abundance} --output {outpath}".split()
    r = runner.invoke(main, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_metacoag_megahit_run(runner, tmp_dir):
    outpath = tmp_dir
    dir_name = DATADIR / "5G_MEGAHIT"
    graph = dir_name / "final.gfa"
    contigs = dir_name / "final.contigs.fa"
    abundance = dir_name / "abundance.tsv"
    args = f"--assembler megahit --graph {graph} --contigs {contigs} --abundance {abundance} --output {outpath}".split()
    r = runner.invoke(main, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_continue_option_is_available(runner):
    result = runner.invoke(main, ["--help"])

    assert result.exit_code == 0
    assert "--continue" in result.output


def test_bidirectional_map_can_be_checkpointed(tmp_path):
    mapping = BidirectionalMap()
    mapping[1] = "contig"
    checkpoint_path = tmp_path / "mapping.pickle"

    with open(checkpoint_path, "wb") as handle:
        pickle.dump(mapping, handle)

    with open(checkpoint_path, "rb") as handle:
        restored = pickle.load(handle)

    assert restored[1] == "contig"
    assert restored.inverse["contig"] == 1


def test_continue_resumes_after_last_completed_stage(tmp_path, monkeypatch):
    input_paths = {}
    for name in ("graph", "contigs", "abundance", "hmm"):
        path = tmp_path / name
        path.write_text(name)
        input_paths[name] = str(path)

    args = SimpleNamespace(
        assembler="custom",
        graph=input_paths["graph"],
        contigs=input_paths["contigs"],
        abundance=input_paths["abundance"],
        paths=None,
        output=str(tmp_path / "output"),
        hmm=input_paths["hmm"],
        prefix="resume",
        min_length=1000,
        p_intra=0.1,
        p_inter=0.01,
        d_limit=20,
        depth=10,
        n_mg=108,
        no_cut_tc=False,
        mg_threshold=0.5,
        bin_mg_threshold=0.33333,
        min_bin_size=200000,
        delimiter=",",
        nthreads=2,
        continue_run=False,
    )
    calls = []

    monkeypatch.setattr(metacoag_pipeline, "log_inputs", lambda *unused: None)
    monkeypatch.setattr(
        metacoag_pipeline,
        "load_graph",
        lambda *unused: (calls.append("graph") or "graph-data", "isolated"),
    )
    monkeypatch.setattr(
        metacoag_pipeline,
        "extract_features",
        lambda *unused: calls.append("features") or "feature-data",
    )
    monkeypatch.setattr(
        metacoag_pipeline,
        "parse_markers",
        lambda *unused: calls.append("markers") or "marker-data",
    )
    monkeypatch.setattr(
        metacoag_pipeline,
        "match_seed_bins",
        lambda *unused: calls.append("seed_bins") or {"seed": True},
    )

    def interrupt_propagation(*unused):
        calls.append("propagated_bins")
        raise RuntimeError("interrupted")

    monkeypatch.setattr(
        metacoag_pipeline,
        "propagate_bins",
        interrupt_propagation,
    )

    with pytest.raises(RuntimeError, match="interrupted"):
        metacoag_runner.run(args)

    config = metacoag_pipeline.PipelineConfig.from_args(args)
    checkpoint_path = metacoag_runner._checkpoint_path(config)
    assert checkpoint_path.is_file()
    assert calls == [
        "graph",
        "features",
        "markers",
        "seed_bins",
        "propagated_bins",
    ]

    calls.clear()
    args.continue_run = True

    def completed_stage_called(*unused):
        raise AssertionError("a completed stage was rerun")

    monkeypatch.setattr(metacoag_pipeline, "load_graph", completed_stage_called)
    monkeypatch.setattr(metacoag_pipeline, "extract_features", completed_stage_called)
    monkeypatch.setattr(metacoag_pipeline, "parse_markers", completed_stage_called)
    monkeypatch.setattr(metacoag_pipeline, "match_seed_bins", completed_stage_called)
    monkeypatch.setattr(
        metacoag_pipeline,
        "propagate_bins",
        lambda *unused: calls.append("propagated_bins") or {"propagated": True},
    )
    monkeypatch.setattr(
        metacoag_pipeline,
        "merge_bins",
        lambda *unused: calls.append("merge_plan") or "merge-plan",
    )
    monkeypatch.setattr(
        metacoag_pipeline,
        "write_output",
        lambda *unused: calls.append("output"),
    )

    metacoag_runner.run(args)

    assert calls == ["propagated_bins", "merge_plan", "output"]
    assert not checkpoint_path.exists()


def test_continue_rejects_changed_options(tmp_path):
    input_paths = {}
    for name in ("graph", "contigs", "abundance", "hmm"):
        path = tmp_path / name
        path.write_text(name)
        input_paths[name] = str(path)

    args = SimpleNamespace(
        assembler="custom",
        graph=input_paths["graph"],
        contigs=input_paths["contigs"],
        abundance=input_paths["abundance"],
        paths=None,
        output=str(tmp_path / "output"),
        hmm=input_paths["hmm"],
        prefix="",
        min_length=1000,
        p_intra=0.1,
        p_inter=0.01,
        d_limit=20,
        depth=10,
        n_mg=108,
        no_cut_tc=False,
        mg_threshold=0.5,
        bin_mg_threshold=0.33333,
        min_bin_size=200000,
        delimiter=",",
        nthreads=2,
        continue_run=False,
    )
    config = metacoag_pipeline.PipelineConfig.from_args(args)
    metacoag_pipeline.validate_inputs(config)
    metacoag_runner._save_checkpoint(config, "graph", {"graph_data": "graph"})

    args.continue_run = True
    args.min_length = 2000
    changed_config = metacoag_pipeline.PipelineConfig.from_args(args)

    with pytest.raises(ValueError, match="does not match"):
        metacoag_runner._load_checkpoint(changed_config)
