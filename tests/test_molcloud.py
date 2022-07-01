import molcloud
import matplotlib.pyplot as plt


def test_plot_molcloud():
    with open("tests/test.smi", "r") as f:
        smls = f.read().splitlines()
    plt.figure(figsize=(10, 3))
    molcloud.plot_molcloud(smls[:1000])
    plt.savefig("test.png")
    plt.close()
