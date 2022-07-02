from .version import __version__
from rdkit.rdBase import BlockLogs
import networkx as nx
import pygraphviz
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
import tqdm

BlockLogs()

_atom_colors = {
    6: "#444444",
    15: "#1BBC9B",
    17: "#a895bb",
    8: "#F06060",
    16: "#F3B562",
    7: "#80cedb",
    35: "#a895bb",
    53: "#a895bb",
}
_unknown_color = "#AAAAAA"
_background_color = "#f5f4e9"


def _custom_layout(G, prog, ratio, args=''):

    A = nx.nx_agraph.to_agraph(G)
    A.graph_attr.update(ratio=ratio)
    # A.write("file.dot")
    A.layout(prog=prog, args=args)
    node_pos = {}
    for n in G:
        node = pygraphviz.Node(A, n)
        try:
            xs = node.attr["pos"].split(",")
            node_pos[n] = tuple(float(x) for x in xs)
        except:
            print("no position for node", n)
            node_pos[n] = (0.0, 0.0)
    return node_pos


def _smiles2graph(sml):
    m = Chem.MolFromSmiles(sml)
    if m is None:
        return None
    # m = rdkit.Chem.AddHs(m)
    G = nx.Graph()
    for a in m.GetAtoms():
        G.add_node(a.GetIdx(), element=a.GetAtomicNum())
    for j in m.GetBonds():
        u = j.GetBeginAtomIdx()
        v = j.GetEndAtomIdx()
        G.add_edge(u, v)
    return G


def _colors(G):
    colors = []
    for n, d in G.nodes(data=True):
        try:
            c = _atom_colors[d["element"]]
        except KeyError as e:
            c = _unknown_color
        colors.append(c)
    return colors


def plot_molcloud(examples, background_color=_background_color, node_size=10, quiet=False):
    G = None
    for smi in tqdm.tqdm(examples, disable=quiet):
        g = _smiles2graph(smi)
        if g is None:
            continue
        if G is None:
            G = g
        else:
            G = nx.disjoint_union(g, G)
    c = _colors(G)
    fig = plt.gcf()
    ratio = fig.get_figheight() / fig.get_figwidth()
    pos = _custom_layout(G, prog="neato", ratio=ratio,
                         args="-Gmaxiter=5000 -Gepsilon=0.00001")
    nx.draw(G, pos, node_size=node_size, node_color=c)
    ax = plt.gca()
    ax.set_facecolor(background_color)
    ax.axis("off")
    fig.set_facecolor(background_color)
