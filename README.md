# molcloud [![build](https://github.com/whitead/molcloud/actions/workflows/tests.yml/badge.svg)](https://whitead.github.io/molcloud/)[![PyPI version](https://badge.fury.io/py/molcloud.svg)](https://badge.fury.io/py/molcloud)

This package draws molecules (or RNA, thanks to [@Shunsuke-1994](https://github.com/Shunsuke-1994)) in a big canvas packed together. See examples below.

## Install

Make sure you have [pygraphviz installed](https://pygraphviz.github.io/documentation/stable/install.html)

```sh
pip install molcloud
```

## Usage

```sh
molcloud [smiles-file]

```

and the output will go to `cloud.png`. You can specify a few things too:

```sh
molcloud [smiles-file] --output-file [output-file] --width 10 --node-size 25
```

To embed the resulting image in some custom shape, use the flag --template
```sh
molcloud [smiles-file] --output-file [output-file] --template [template-file]
```

Use `molcloud --help` for complete options. `smiles-file` should contain smiles, one per line like:

```plain
O=C(OC)C=1C=CC2=NC=C(C(=O)OCC)C(NCC(O)C)=C2C1
O=C1C2=NC=CC3=C(OC)C=4OCOC4C(C=5C=C(OC)C(OC)=C(OC)C15)=C23
```

Adjust width as you add more molecules. The drawing is always square (sorry).

## RNA Install

Thanks to [@Shunsuke-1994](https://github.com/Shunsuke-1994)! To install layout RNA, install the extra packages:

```sh
pip install molcloud[all]
```

## RNA Usage

```sh
rnacloud [fasta-file]
```

where `fasta-file` should contain sequence and bracket notations, three lines per 1 sequence like:
```
>seq_0
UUCCAGCACCUGAUGUUCGAAUUUAAAUCGGCUCAACGAG
(((.((((.....)))).)))......(((......))).
```

## Molecule Example
![test](https://user-images.githubusercontent.com/908389/176980703-bc814295-ee37-4c41-a31b-6b75bb420659.png)

## Example with template

```sh
molcloud tests/test.smi --template tests/shapes/beaker.png --width 15
```

## RNA Example

![rna](https://user-images.githubusercontent.com/908389/177061306-8caea628-12a4-4ccd-ae7d-ae240ba3adb1.png)


## Animation Example
