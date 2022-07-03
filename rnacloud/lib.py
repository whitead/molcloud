from .version import __version__
import networkx as nx
import pygraphviz
import matplotlib.pyplot as plt
import tqdm
import forgi.graph.bulge_graph as fgb

_atom_colors = {
    "A": "#61D04F",
    "C": "#2297E6",
    "G": "#F5C710",
    "U": "#DF536B"}

_bond_colors = {
    "covalent": "#000000",
    "hydrogen": "#999999"
    }

_unknown_color = "#AAAAAA"
_background_color = "#f5f4e9"

def _fasta_text2graph(fasta_text):
    bg = fgb.BulgeGraph.from_fasta_text(fasta_text)[0]
    G = nx.Graph()
    residues = []
    for d in bg.defines:
        prev = None

        for r in bg.define_residue_num_iterator(d):
            G.add_node(r, nucleotide = str(bg._seq)[r-1])
            residues += [r]

    # Add links along the backbone
    residues.sort()
    prev = None
    for r in residues:
        if prev is not None:
            G.add_edge(prev, r, bond = "covalent")
        prev = r

    # Add links along basepairs
    for s in bg.stem_iterator():
        for (f, t) in bg.stem_bp_iterator(s):
            G.add_edge(f, t, bond = "hydrogen")

    return G


def _colors(G):
    node_colors = []
    for n, d in G.nodes(data=True):
        try:
            c = _atom_colors[d["nucleotide"]]
        except KeyError as e:
            c = _unknown_color
        node_colors.append(c)
        
    edge_colors = []
    for *e, d in G.edges(data=True):
        try:
            c = _bond_colors[d["bond"]]
        except KeyError as e:
            c = _unknown_color
        edge_colors.append(c)
    return node_colors, edge_colors 


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


def plot_rnacloud(fasta_texts, background_color=background_color, node_size=10, quiet=False):
    G = None

    for fasta_text in tqdm.tqdm(fasta_texts, disable=quiet):
        g = _fasta_text2graph(fasta_text)
        if g is None:
            continue
        if G is None:
            G = g
        else:
            G = nx.disjoint_union(g, G)
    node_colors, edge_colors = _colors(G)
    fig = plt.gcf()
    ratio = fig.get_figheight() / fig.get_figwidth()
    pos = _custom_layout(G, prog="neato", ratio=ratio,
                         args="-Gmaxiter=5000 -Gepsilon=0.00001")
    nx.draw(G, pos, node_size=node_size, node_color=node_colors, edge_color = edge_colors)
    ax = plt.gca()
    ax.set_facecolor(background_color)
    ax.axis("off")
    fig.set_facecolor(background_color)