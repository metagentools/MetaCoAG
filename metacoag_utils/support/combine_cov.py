#!/usr/bin/python3

"""combine_cov.py: Combine multiple coverage files of samples from CoverM.
"""

import click
import glob

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
    """
    combine_cov: Combine multiple coverage files of samples from CoverM
    """

    # Get coverage values from samples
    # ---------------------------------------------------

    # Get coverage files
    cov_files = glob.glob(f"{covpath}/*.tsv")

    final_df = pd.DataFrame()

    for file in cov_files:
        df = pd.read_csv(file, sep="\t", header=0)

        if final_df.empty:
            final_df = df
        else:
            final_df = pd.concat(
                [final_df, df[list(df.columns)[1]]], axis=1, join="inner"
            )

    print(f"Dataframe shape: {final_df.shape}")

    # Save dataframe to file
    final_df.to_csv(f"{output}coverage.tsv", sep="\t", index=False, header=False)
    final_df.to_csv(
        f"{output}coverage_with_header.tsv", sep="\t", index=False, header=True
    )
    print(f"The combined coverage values can be found at {output}coverage.tsv")

    # Exit program
    # --------------

    print("Thank you for using combine_cov!")


if __name__ == "__main__":
    main()
