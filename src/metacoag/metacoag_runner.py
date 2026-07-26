#!/usr/bin/env python3

import dataclasses
import logging
import os
import pathlib
import pickle
import tempfile
import time

from metacoag import metacoag_pipeline


__author__ = "Vijini Mallawaarachchi and Yu Lin"
__copyright__ = "Copyright 2020, MetaCoAG Project"
__license__ = "GPL-3.0"
__version__ = "1.3.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Stable Release"

CHECKPOINT_VERSION = 1
STAGE_ORDER = {
    "graph": 1,
    "features": 2,
    "markers": 3,
    "seed_bins": 4,
    "propagated_bins": 5,
    "merge_plan": 6,
}


def _configure_logger(config):
    logger = logging.getLogger(f"MetaCoAG {__version__}")
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    console_handler.setLevel(logging.INFO)

    log_path = config.output / f"{config.prefix}metacoag.log"
    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    return logger, (console_handler, file_handler)


def _checkpoint_path(config):
    return config.output / f".{config.prefix}metacoag_checkpoint.pkl"


def _file_signature(path):
    if path is None:
        return None

    file_path = pathlib.Path(path)
    if not file_path.is_file():
        return {"path": str(file_path.resolve()), "missing": True}

    stat = file_path.stat()
    return {
        "path": str(file_path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _config_signature(config):
    values = dataclasses.asdict(config)
    values.pop("continue_run")
    values.pop("nthreads")
    values["output"] = str(config.output.resolve())

    for name in ("graph", "contigs", "abundance", "paths", "hmm"):
        values[name] = _file_signature(getattr(config, name))

    return values


def _save_checkpoint(config, stage, payload):
    checkpoint_path = _checkpoint_path(config)
    checkpoint = {
        "version": CHECKPOINT_VERSION,
        "stage": stage,
        "signature": _config_signature(config),
        "payload": payload,
    }

    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=config.output,
            prefix=f"{checkpoint_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = pathlib.Path(handle.name)
            pickle.dump(checkpoint, handle, protocol=pickle.HIGHEST_PROTOCOL)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, checkpoint_path)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def _load_checkpoint(config):
    checkpoint_path = _checkpoint_path(config)
    if not checkpoint_path.is_file():
        raise FileNotFoundError(
            f"No MetaCoAG checkpoint found in {config.output}. "
            "Run without --continue to start a new analysis."
        )

    try:
        with open(checkpoint_path, "rb") as handle:
            checkpoint = pickle.load(handle)
    except Exception as error:
        raise RuntimeError(
            f"Could not read MetaCoAG checkpoint: {checkpoint_path}"
        ) from error

    if (
        not isinstance(checkpoint, dict)
        or checkpoint.get("version") != CHECKPOINT_VERSION
        or checkpoint.get("stage") not in STAGE_ORDER
        or not isinstance(checkpoint.get("payload"), dict)
    ):
        raise RuntimeError(
            f"Unsupported or corrupt MetaCoAG checkpoint: {checkpoint_path}"
        )

    if checkpoint.get("signature") != _config_signature(config):
        raise ValueError(
            "The checkpoint does not match the current inputs or options. "
            "Run without --continue to start a new analysis."
        )

    return checkpoint


def _stage_completed(checkpoint, stage):
    return (
        checkpoint is not None
        and STAGE_ORDER[checkpoint["stage"]] >= STAGE_ORDER[stage]
    )


def run(args):
    config = metacoag_pipeline.PipelineConfig.from_args(args)
    metacoag_pipeline.validate_inputs(config)
    logger, handlers = _configure_logger(config)
    start_time = time.time()

    try:
        metacoag_pipeline.log_inputs(config, logger)
        logger.info("MetaCoAG started")

        checkpoint_path = _checkpoint_path(config)
        if config.continue_run:
            checkpoint = _load_checkpoint(config)
            payload = checkpoint["payload"]
            logger.info(
                "Continuing from checkpoint after stage: %s",
                checkpoint["stage"],
            )
        else:
            checkpoint_path.unlink(missing_ok=True)
            checkpoint = None
            payload = {}

        if _stage_completed(checkpoint, "graph"):
            graph_data = payload["graph_data"]
            isolated = payload["isolated"]
        else:
            graph_data, isolated = metacoag_pipeline.load_graph(config, logger)
            payload.update(graph_data=graph_data, isolated=isolated)
            _save_checkpoint(config, "graph", payload)

        if _stage_completed(checkpoint, "features"):
            features = payload["features"]
        else:
            features = metacoag_pipeline.extract_features(
                config, graph_data, isolated, logger
            )
            payload["features"] = features
            _save_checkpoint(config, "features", payload)

        if _stage_completed(checkpoint, "markers"):
            markers = payload["markers"]
        else:
            markers = metacoag_pipeline.parse_markers(
                config, graph_data, features, logger
            )
            payload["markers"] = markers
            _save_checkpoint(config, "markers", payload)

        if _stage_completed(checkpoint, "seed_bins"):
            bins = payload["bins"]
        else:
            bins = metacoag_pipeline.match_seed_bins(
                config, graph_data, features, markers, logger
            )
            payload["bins"] = bins
            _save_checkpoint(config, "seed_bins", payload)

        if not _stage_completed(checkpoint, "propagated_bins"):
            bins = metacoag_pipeline.propagate_bins(
                config, graph_data, features, markers, bins, logger
            )
            payload["bins"] = bins
            _save_checkpoint(config, "propagated_bins", payload)

        if _stage_completed(checkpoint, "merge_plan"):
            merge_plan = payload["merge_plan"]
        else:
            merge_plan = metacoag_pipeline.merge_bins(config, features, bins, logger)
            payload["merge_plan"] = merge_plan
            _save_checkpoint(config, "merge_plan", payload)

        metacoag_pipeline.write_output(config, graph_data, bins, merge_plan, logger)
        checkpoint_path.unlink(missing_ok=True)

        logger.info("Elapsed time: %s seconds", time.time() - start_time)
        logger.info("Thank you for using MetaCoAG!")
    finally:
        for handler in handlers:
            logger.removeHandler(handler)
            handler.close()


def main(args):
    run(args)


if __name__ == "__main__":
    main()
