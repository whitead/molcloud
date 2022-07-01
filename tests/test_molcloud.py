import molcloud
import matplotlib.pyplot as plt


def test_plot_molcloud():
    with open("tests/test.smi", "r") as f:
        smls = f.read().splitlines()
    plt.figure(figsize=(15, 15))
    molcloud.plot_molcloud(smls[:500])
    plt.savefig("test.png")
    plt.close()
