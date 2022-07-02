from importlib.resources import path
from .lib import plot_molcloud
import matplotlib.pyplot as plt
import click


@click.command()
@click.argument("smiles-file", type=click.Path(exists=True))
@click.option("--output-file", default="cloud.png", type=click.Path(exists=False), help="Path to the output file")
@click.option("--width", type=int, default=10, help="Width of the output image")
@click.option("--background-color", default="#f5f4e9", help="Background color of the output image")
@click.option("--node-size", type=int, default=10, help="Size of the nodes")
@click.option("--quiet", is_flag=True, default=False, help="Don't show progress")
def main(smiles_file, output_file, width, background_color, node_size, quiet):
    with open(smiles_file, "r") as f:
        smls = f.read().splitlines()
    plt.figure(figsize=(width, width))
    plot_molcloud(smls, background_color=background_color,
                  node_size=node_size, quiet=quiet)
    plt.savefig(output_file)
