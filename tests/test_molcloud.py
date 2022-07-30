import molcloud
import random
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def dtest_custom_layout():
    # make random graph
    G = nx.gnp_random_graph(100, 0.5)
    pos = molcloud.lib.custom_layout(G, 'neato')
    new_pos = molcloud.lib.custom_layout(G, 'neato', start_pos=pos)
    # make sure it didn't change
    # doesn't always hold, because of graphviz
    # for n in G:
    #    np.testing.assert_allclose(new_pos[n], pos[n], atol=0.1)


def dtest_animate_molcloud():
    with open("tests/test.smi", "r") as f:
        smls = f.read().splitlines()
    plt.figure(figsize=(7, 7))
    molcloud.darkmode()
    molcloud.animate_molcloud(smls[:15], prog='sfdp', background_color='#333')


def test_plot_molcloud():
    with open("tests/test.smi", "r") as f:
        smls = f.read().splitlines()
    plt.figure(figsize=(20, 20))
    molcloud.plot_molcloud(smls[:2])
    plt.savefig("cover.png")
    plt.close()


def test_plot_molcloud_template():
    with open("tests/test.smi", "r") as f:
        smls = f.read().splitlines()
    plt.figure(figsize=(15, 15))
    molcloud.plot_molcloud(
        smls, template="tests/shapes/beaker.png", thresh=0.2)
    plt.savefig("beaker.png")
    plt.close()


def test_plot_rnacloud():
    fasta_texts = []
    fasta_text = ""
    with open('tests/test.fa') as f:
        for line in f.readlines():
            if line.startswith(">"):
                if fasta_text != "":
                    fasta_texts.append(fasta_text)
                fasta_text = ""
            else:
                fasta_text += line
    plt.figure(figsize=(15, 15))
    molcloud.plot_rnacloud(fasta_texts)
    plt.savefig("rna.png")
    plt.close()
