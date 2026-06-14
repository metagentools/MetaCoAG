#!/usr/bin/python3

"""Combine multiple CoverM sample coverage files."""

from pathlib import Path

import click
import pandas as pd


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2020, MetaCoAG Project"
__license__ = "GPL-3.0"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"


@click.command()
@click.option(
    "--covpath",
    help="path to the .tsv files from CoverM",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--output",
    help="path to the output folder",
    type=click.Path(dir_okay=True, writable=True, readable=True),
    required=True,
)
def main(covpath, output):
    """Combine multiple coverage files from CoverM."""
    coverage_files = sorted(Path(covpath).glob("*.tsv"))
    combined_coverage = pd.DataFrame()

    for coverage_file in coverage_files:
        coverage_frame = pd.read_csv(coverage_file, sep="\t", header=0)

        if combined_coverage.empty:
            combined_coverage = coverage_frame
        else:
            sample_column = coverage_frame.columns[1]
            combined_coverage = pd.concat(
                [combined_coverage, coverage_frame[sample_column]],
                axis=1,
                join="inner",
            )

    print(f"Dataframe shape: {combined_coverage.shape}")

    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)
    coverage_path = output_path / "coverage.tsv"
    header_path = output_path / "coverage_with_header.tsv"

    combined_coverage.to_csv(coverage_path, sep="\t", index=False, header=False)
    combined_coverage.to_csv(header_path, sep="\t", index=False, header=True)
    print(f"The combined coverage values can be found at {coverage_path}")
    print("Thank you for using combine_cov!")


if __name__ == "__main__":
    main()
