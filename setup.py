import os
from glob import glob
from typing import Optional
from setuptools import setup

exec(open("molcloud/version.py").read())

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="molcloud",
    version=__version__,
    description="Insert SVGs into matplotlib figures",
    author="Andrew White",
    author_email="andrew.white@rochester.edu",
    url="https://github.com/whitead/molcloud",
    license="MIT",
    packages=["molcloud"],
    install_requires=["networkx", "matplotlib",
                      "pygraphviz", "rdkit", "click", "tqdm"],
    optional_dependencies={"rna": ["forgi"]},
    test_suite="tests",
        entry_points="""
        [console_scripts]
        molcloud=molcloud.main:smiles
        rnacloud=molcloud.main:rna
        """,
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Framework :: Matplotlib",
    ],
)
