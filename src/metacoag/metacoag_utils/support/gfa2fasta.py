"""Convert Flye and Miniasm GFA segment sequences to FASTA.

The Flye assembly graph file (assembly_graph.gfa) is expected as input.
"""

import logging
import re

from pathlib import Path

import click

from cogent3.format.fasta import alignment_to_fasta


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2020, MetaCoAG Project"
__license__ = "GPL-3.0"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"


@click.command()
@click.option(
    "--graph",
    help="path to the assembly graph file",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--output",
    help="path to the output folder",
    type=click.Path(dir_okay=True, writable=True, readable=True),
    required=True,
)
@click.option(
    "--log",
    help="path to log file",
    type=click.Path(dir_okay=True, writable=True, readable=True),
    required=False,
)
def main(graph, output, log):
    assembly_graph_path = Path(graph)
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("metacoag.gfa2fasta")
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    console_handler.setLevel(logging.INFO)
    logger.addHandler(console_handler)

    log_path = Path(log) if log is not None else output_path / "gfa2fasta.log"
    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.info("Obtaining edge sequences")

    sequences = {}
    with open(assembly_graph_path) as graph_file:
        for line in graph_file:
            if not line.startswith("S"):
                continue
            fields = line.split("\t")
            sequences[fields[1]] = re.sub("[^GATC]", "", fields[2].upper())

    logger.info("Writing edge sequences to FASTA file")

    fasta_path = output_path / "edges.fasta"
    with open(fasta_path, "w") as output_file:
        output_file.write(alignment_to_fasta(sequences))

    logger.info("The FASTA file with unitig sequences can be found at %s", fasta_path)
    logger.info("Thank you for using gfa2fasta!")

    logger.removeHandler(file_handler)
    logger.removeHandler(console_handler)
    file_handler.close()
    console_handler.close()


if __name__ == "__main__":
    main()
