"""gfa2fasta.py: Obtain the sequences corresponding to edges in the Flye and Miniasm assembly graphs in FASTA format.

The assembly graph file of Flye (assembly_graph.gfa) should be provided as inputs.

"""

import click
import logging
import os
import pathlib
import re
import sys

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
    # Get arguments
    # -----------------------

    assembly_graph_file = graph
    output_path = pathlib.Path(output)
    output_path.mkdir(parents=True, exist_ok=True)
    log_file = log
    prefix = ""

    # Setup logger
    # ----------------------------------------------------------------------

    logger = logging.getLogger("gfa2fasta")
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    # Setup output path for log file
    if log_file is None:
        fileHandler = logging.FileHandler(output_path / "gfa2fasta.log")
    else:
        fileHandler = logging.FileHandler(f"{log_file}")

    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    # Check assembly graph file
    if not os.path.isfile(assembly_graph_file):
        logger.error(
            "Failed to open the assembly graph file. Please make sure to provife the .gfa file."
        )
        logger.info("Exiting gfa2fasta.py...\nBye...!\n")
        sys.exit(1)

    # Get the sequences corresponding to edges of the graph.
    # ---------------------------------------------------

    logger.info("Obtaining edge sequences")

    seqs = {}

    with open(assembly_graph_file) as file:
        line = file.readline()

        while line != "":
            if "S" in line:
                strings = line.split("\t")
                seqs[str(strings[1])] = re.sub("[^GATC]", "", str(strings[2]).upper())

            line = file.readline()

    logger.info("Writing edge sequences to FASTA file")

    with open(output_path / f"{prefix}edges.fasta", "w") as output:
        output.write(alignment_to_fasta(seqs))

    logger.info(
        f"The FASTA file with unitig sequences can be found at {output.name}"
    )

    # Exit program
    # --------------

    logger.info("Thank you for using gfa2fasta!")


if __name__ == "__main__":
    main()
