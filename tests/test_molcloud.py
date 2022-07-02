import molcloud
import random
import matplotlib.pyplot as plt


def test_plot_molcloud():
    with open("tests/test.smi", "r") as f:
        smls = f.read().splitlines()
    random.shuffle(smls)
    plt.figure(figsize=(20, 20))
    molcloud.plot_molcloud(smls)
    plt.savefig("cover.png")
    plt.close()
