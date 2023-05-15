#!/usr/bin/env python3

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

packages = setuptools.find_packages()
package_data = {"metacoag_utils": ["metacoag_utils/*", "metacoag_utils/auxiliary/*"]}

data_files = [(".", ["LICENSE", "README.md"])]

setuptools.setup(
    name="metacoag",
    version="1.1.2",
    zip_safe=True,
    author="Vijini Mallawaarachchi and Yu Lin",
    author_email="viji.mallawaarachchi@gmail.com",
    description="MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/metagentools/MetaCoAG",
    license="GPL-3.0",
    packages=packages,
    package_data=package_data,
    data_files=data_files,
    include_package_data=True,
    scripts=["metacoag"],
    entry_points={
        "console_scripts": [
            "combine_cov=metacoag_utils.support.combine_cov:main",
        ],
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "biopython",
        "python-igraph",
        "networkx",
        "scipy",
        "numpy",
        "tqdm",
        "pandas",
        "hmmer",
        "click"
    ],
    python_requires=">=3.7",
)
