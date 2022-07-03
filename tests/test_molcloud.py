import molcloud
import random
import matplotlib.pyplot as plt


def test_plot_molcloud():
    with open("tests/test.smi", "r") as f:
        smls = f.read().splitlines()
    random.shuffle(smls)
    plt.figure(figsize=(7, 7))
    molcloud.plot_molcloud(smls[:25])
    plt.savefig("cover.png")
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
