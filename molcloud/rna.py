import forgi.graph.bulge_graph as fgb
from .lib import *

_rna_colors = {
    "A": "#61D04F",
    "C": "#2297E6",
    "G": "#F5C710",
    "U": "#DF536B"}

_rna_bond_colors = {
    "covalent": "#000000",
    "hydrogen": "#999999"
}


def _fasta_text2graph(fasta_text):
    bg = fgb.BulgeGraph.from_fasta_text(fasta_text)[0]
    G = nx.Graph()
    residues = []
    for d in bg.defines:
        prev = None

        for r in bg.define_residue_num_iterator(d):
            G.add_node(r, nucleotide=str(bg._seq)[r-1])
            residues += [r]

    # Add links along the backbone
    residues.sort()
    prev = None
    for r in residues:
        if prev is not None:
            G.add_edge(prev, r, bond="covalent")
        prev = r

    # Add links along basepairs
    for s in bg.stem_iterator():
        for (f, t) in bg.stem_bp_iterator(s):
            G.add_edge(f, t, bond="hydrogen")

    return G


def _colors(G):
    node_colors = []
    for n, d in G.nodes(data=True):
        try:
            c = _rna_colors[d["nucleotide"]]
        except KeyError as e:
            c = unknown_color
        node_colors.append(c)

    edge_colors = []
    for *e, d in G.edges(data=True):
        try:
            c = _rna_bond_colors[d["bond"]]
        except KeyError as e:
            c = unknown_color
        edge_colors.append(c)
    return node_colors, edge_colors


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
    plot_graphs(G, node_colors, edge_colors, background_color, node_size)
