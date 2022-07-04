from .version import __version__
from rdkit.rdBase import BlockLogs
import networkx as nx
import pygraphviz
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
import tqdm
import cv2
import pandas as pd

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


def _smiles2graph(sml, aim=False, molid=0, offset=0):

    # Track atoms in molecule
    df = []
    m = Chem.MolFromSmiles(sml)
    if m is None:
        return None
    # m = rdkit.Chem.AddHs(m)
    G = nx.Graph()
    for i,a in enumerate(m.GetAtoms()):
        idx=a.GetIdx()
        G.add_node(idx, element=a.GetAtomicNum())
        if aim:    df.append([molid,idx+offset,1])
    for j in m.GetBonds():
        u = j.GetBeginAtomIdx()
        v = j.GetEndAtomIdx()
        G.add_edge(u, v)

    if aim:
        df = pd.DataFrame(df, columns=["molid","atom","isin"])
        return G,df
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

def _loadTemplate(template):

    ########    Do some image processing on template image.
    img = cv2.imread(template,0)
    # Canny edge detection
    img = cv2.Canny(img,100,200)
    ctrs = cv2.findContours(img, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)

    # Get largest contour
    max_c=[]
    for c in ctrs[0]:
        if len(c)>len(max_c):
            max_c=c

    M = cv2.moments(max_c)
    x = int(M["m10"] / M["m00"])
    y = int(M["m01"] / M["m00"])

    # Initiale mask for flood filling
    width, height = img.shape
    mask = img2 = np.ones((width + 2, height + 2), np.uint8) * 255
    mask[1:width, 1:height] = 0

    # Generate intermediate image, draw largest contour, flood filled
    temp = np.zeros(img.shape, np.uint8)
    temp = cv2.drawContours(temp, max_c, -1, 1, cv2.FILLED)
    _, temp, mask, _ = cv2.floodFill(temp, mask, (x, y), 1)
    return temp

def _dropMols(G, pos, c, temp, moldf):

    # Get pos in same scale as template
    posdf = pd.DataFrame(pos)
    posdf = posdf.divide(posdf.max(axis=1),axis=0)
    posdf.iloc[1] *= temp.shape[0] - 1
    posdf.iloc[0] *= temp.shape[1] - 1
    posdf = posdf.astype(int).to_dict('list')

    # First: label atoms as in or out boundary
    for i,atom in enumerate(G):
        coords = posdf[atom]
        if temp[temp.shape[0]-coords[1]-1, coords[0]] == 0:
            moldf.loc[i, "isin"] = 0

    # Decide what molecules to drop
    gb = moldf.groupby("molid").apply(lambda x: x.sum()/x.count())
    gb = (gb["isin"]>0.7).drop(columns="atom")

    moldf = pd.merge(moldf.drop(columns="isin"), gb, left_on="molid", right_index=True)

    new_c = []
    print(moldf.head(50))
    for at in G:
        print(at)

    for i in moldf.index:
        isin, atom = moldf.loc[i, ["isin","atom"]]
        if not isin:
            G.remove_node(atom)
            pos.pop(atom)
        else:
            new_c.append(c[atom])
            
    return G, pos, new_c


def plot_molcloud(examples, background_color=_background_color, node_size=10, quiet=False,
                  template=None, repeat=0):


    # repeat dataset N times, so image is more filled.
    # But reshuffle so locations are not so close
    for _ in range(repeat):
        shuf_examp = examples.copy()
        np.random.shuffle(shuf_examp)
        examples += list(shuf_examp)


    G = None

    # Store for each atom in each molecule: is in or out the region.
    moldf = pd.DataFrame(columns=["molid","atom","isin"])

    for i, smi in enumerate(tqdm.tqdm(examples, disable=quiet)):
        if template:
            g, df_ = _smiles2graph(smi, aim=True, molid=i, offset=moldf.shape[0])
            moldf = pd.concat([moldf, df_])
        else:
            g = _smiles2graph(smi)

        if g is None:
            continue
        if G is None:
            G = g
        else:
            G = nx.disjoint_union(G, g)


    moldf = moldf.reset_index(drop=True)
    fig = plt.gcf()

   
    if template:    # Get a binary classified map
        temp = _loadTemplate(template)
        ratio = temp.shape[0] / temp.shape[1]
    else:
        ratio = fig.get_figheight() / fig.get_figwidth()

    c = _colors(G)
    pos = _custom_layout(G, prog="neato", ratio=ratio,
                         args="-Gmaxiter=5000 -Gepsilon=0.00001")

    if template:
        ### Remove molecules outside the white region
        G, pos, c = _dropMols(G, pos, c, temp, moldf)
        


    nx.draw(G, pos, node_size=node_size, node_color=c)
    ax = plt.gca()
    ax.set_facecolor(background_color)
    ax.axis("off")
    fig.set_facecolor(background_color)
