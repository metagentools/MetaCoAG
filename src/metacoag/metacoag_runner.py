#!/usr/bin/env python3

import logging
import time

from metacoag import metacoag_pipeline

__author__ = "Vijini Mallawaarachchi and Yu Lin"
__copyright__ = "Copyright 2020, MetaCoAG Project"
__license__ = "GPL-3.0"
__version__ = "1.2.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Stable Release"


def _configure_logger(config):
    logger = logging.getLogger(f"MetaCoaAG {__version__}")
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


def run(args):
    config = metacoag_pipeline.PipelineConfig.from_args(args)
    metacoag_pipeline.validate_inputs(config)
    logger, handlers = _configure_logger(config)
    start_time = time.time()

    try:
        metacoag_pipeline.log_inputs(config, logger)
        logger.info("MetaCoAG started")

        graph_data, isolated = metacoag_pipeline.load_graph(config, logger)
        features = metacoag_pipeline.extract_features(
            config, graph_data, isolated, logger
        )
        markers = metacoag_pipeline.parse_markers(
            config, graph_data, features, logger
        )
        bins = metacoag_pipeline.match_seed_bins(
            config, graph_data, features, markers, logger
        )
        bins = metacoag_pipeline.propagate_bins(
            config, graph_data, features, markers, bins, logger
        )
        merge_plan = metacoag_pipeline.merge_bins(
            config, features, bins, logger
        )
        metacoag_pipeline.write_output(
            config, graph_data, bins, merge_plan, logger
        )

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
