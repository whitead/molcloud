import molcloud
import random
import matplotlib.pyplot as plt


def test_plot_molcloud():
    with open("tests/test.smi", "r") as f:
        smls = f.read().splitlines()
    random.shuffle(smls)
    plt.figure(figsize=(10, 10))
    molcloud.plot_molcloud(smls[:100])
    plt.savefig("cover.png")
    plt.close()
