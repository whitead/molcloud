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
background_color = "#f5f4e9"


def darkmode():
    global background_color
    _atom_colors[6] = '#ef9af2'


def _smooth(data, window_size):
    # smooth with moving average
    # assumes time is 0th axis
    def do_smooth(x):
        cumsum_vec = np.cumsum(np.insert(x, 0, 0))
        ma_vec = (cumsum_vec[window_size:] -
                  cumsum_vec[:-window_size]) / window_size
        return ma_vec
    flat_data = data.reshape(data.shape[0], -1)
    s = np.apply_along_axis(do_smooth, 0, flat_data)
    return s.reshape(-1, *data.shape[1:])


def custom_layout(G, prog, args='', start_pos=None):
    A = nx.nx_agraph.to_agraph(G)
    inputscale = 72  # const in graphviz
    # set start position
    if start_pos is not None:
        for n in G:
            node = pygraphviz.Node(A, n)
            if n in start_pos:
                node.attr["pos"] = f'{start_pos[n][0]:.2f},{start_pos[n][1]:.2f}'
                node.attr["pin"] = True
    A.graph_attr["notranslate"] = True
    A.graph_attr["inputscale"] = inputscale
    A.write('start' + ".dot")
    A.layout(prog=prog, args=args)
    A.write('end' + ".dot")
    node_pos = {}
    for n in G:
        node = pygraphviz.Node(A, n)
        try:
            xs = node.attr["pos"].split(",")
            node_pos[n] = tuple(float(x) for x in xs)
        except:
            print("no position for node", n)
            node_pos[n] = (0.0, 0.0)
    # print maximum node position
    return node_pos


def anim_custom_layout(G, prog, steps=100, iters=3, start_pos=None, window_size=5):
    steps = max(steps, 1)
    do_steps = steps + window_size - 2
    inputscale = 72  # const in graphviz
    A = nx.nx_agraph.to_agraph(G)
    # set start position
    if start_pos is not None:
        for n in G:
            node = pygraphviz.Node(A, n)
            if n in start_pos:
                node.attr["pos"] = f'{start_pos[n][0]:.2f},{start_pos[n][1]:.2f}'
                node.attr["pin"] = False
    A.graph_attr["inputscale"] = inputscale
    A.graph_attr["notranslate"] = True
    node_pos = np.empty((do_steps, len(G), 2))
    for i in range(do_steps):
        A.layout(
            prog=prog, args=f"-Gmaxiter={iters} -Gepsilon=0.00001")
        for j, n in enumerate(G):
            node = pygraphviz.Node(A, n)
            try:
                xs = node.attr["pos"].split(",")
                node_pos[i, j] = tuple(float(x) for x in xs)
            except:
                print("no position for node", n)
                node_pos[i, j] = (0.0, 0.0)
    # now convert to list of dicts
    snode_pos = _smooth(node_pos, window_size)
    pos = []
    for i in range(steps - 1):
        pos.append(dict(zip(G, snode_pos[i])))
    # want final one to be real, not smoothed
    pos.append(dict(zip(G, node_pos[-1])))
    return pos


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


def plot_graphs(G, node_colors, edge_colors, background_color, node_size):
    fig = plt.gcf()
    # get fig width
    #size = fig.get_figwidth(), fig.get_figheight()
    pos = custom_layout(G, prog="neato",
                        args="-Gmaxiter=5000 -Gepsilon=0.00001")
    nx.draw(G, pos, node_size=node_size,
            node_color=node_colors, edge_color=edge_colors)
    ax = plt.gca()
    ax.set_facecolor(background_color)
    ax.axis("off")
    fig.set_facecolor(background_color)


def plot_molcloud(examples, background_color=background_color, node_size=10, quiet=False):
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
    plot_graphs(G, c, '#333', background_color, node_size)


def animate_molcloud(examples, prog="neato", background_color=background_color, node_size=10, quiet=False, duration=3, fps=20):
    from moviepy.editor import VideoClip
    from moviepy.video.io.bindings import mplfig_to_npimage

    all_G = []
    all_c = []
    G = None
    pos = None
    i_steps = 30
    N = 0
    for smi in tqdm.tqdm(examples, disable=quiet):
        g = _smiles2graph(smi)
        if g is None:
            continue
        if G is None:
            G = g
            pos = [anim_custom_layout(G, prog, steps=i_steps)]
        else:
            G = nx.disjoint_union(G, g)
            pos.append(anim_custom_layout(
                G, prog, steps=i_steps, start_pos=pos[-1][-1]))
        all_G.append(G)
        all_c.append(_colors(G))
        N += i_steps
    fig = plt.gcf()

    nodes = nx.draw_networkx_nodes(
        all_G[0], pos[0][0], node_color=all_c[0], node_size=node_size)
    edges = nx.draw_networkx_edges(all_G[0], pos[0][0], edge_color='#FFF')
    ax = plt.gca()
    ax.set_facecolor(background_color)
    ax.axis("off")
    fig.set_facecolor(background_color)
    fig.tight_layout()

    i = 0
    j = 0

    def make_frame(t):
        nonlocal i, j
        # update node positions
        ax.clear()
        nodes = nx.draw_networkx_nodes(
            all_G[i], pos[i][j], node_size=node_size, node_color=all_c[i])
        edges = nx.draw_networkx_edges(all_G[i], pos[i][j], edge_color='#FFF')
        j += 1
        if j >= i_steps:
            j = 0
            i += 1
        i = min(i, len(all_G) - 1)
        return mplfig_to_npimage(fig)

    animation = VideoClip(make_frame, duration=N / fps)
    animation.write_gif('animation.gif', fps=fps)
