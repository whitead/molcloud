# molcloud [![build](https://github.com/whitead/molcloud/actions/workflows/tests.yml/badge.svg)](https://whitead.github.io/molcloud/)[![PyPI version](https://badge.fury.io/py/molcloud.svg)](https://badge.fury.io/py/molcloud)

## Install

```sh
pip install molcloud
```

## Usage

```sh
molcloud [smiles_file] [output_file] --width 10
```

Use `molcloud --help` for complete options. `smiles_file` should contain smiles, one per line like:

```plain
=C(OC)C=1C=CC2=NC=C(C(=O)OCC)C(NCC(O)C)=C2C1
O=C1C2=NC=CC3=C(OC)C=4OCOC4C(C=5C=C(OC)C(OC)=C(OC)C15)=C23
```

Adjust width as you add more molecules. The drawing is always square (sorry).