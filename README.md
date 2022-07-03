# molcloud [![build](https://github.com/whitead/molcloud/actions/workflows/tests.yml/badge.svg)](https://whitead.github.io/molcloud/)[![PyPI version](https://badge.fury.io/py/molcloud.svg)](https://badge.fury.io/py/molcloud)

## Install

Make sure you have [pygraphviz installed](https://pygraphviz.github.io/documentation/stable/install.html)
Also [opencv](https://pypi.org/project/opencv-python/)

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

## Example

![test](https://user-images.githubusercontent.com/908389/176980703-bc814295-ee37-4c41-a31b-6b75bb420659.png)

## Example with template
![test_templ](tests/beaker.png)

