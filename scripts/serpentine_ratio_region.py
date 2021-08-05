# This is a command line utility to plot the serpentine ratio of a region between two input cooler files
# cmdoret, 20210420

import csv
from os.path import basename
from typing import Optional, Tuple
import cooler
import click
import numpy as np
import pandas as pd
import serpentine as serp
import matplotlib.pyplot as plt


def get_ratio(
    clr1: cooler.Cooler,
    clr2: cooler.Cooler,
    region: str,
    no_serp: bool = False,
) -> np.ndarray:
    """Compute the ratio between two input cooler objects in predefined region, use
    serpentine's dynamic binning algorithm and return the resulting log ratio matrix"""
    mat1 = clr1.matrix(balance=False, sparse=False).fetch(region)
    mat2 = clr2.matrix(balance=False, sparse=False).fetch(region)
    # Ratio will be log2(mat1 / mat2)
    if no_serp:
        ratio = np.log2(mat1 / mat2)
        ratio[np.isneginf(ratio)] = np.nan
        ratio[np.isposinf(ratio)] = np.nan
    else:
        _, _, ratio = serp.serpentin_binning(
            mat2, mat1, force_symmetric=True, parallel=12, iterations=50
        )
    return ratio


def get_matrices(
        clr1: cooler.Cooler, clr2: cooler.Cooler, region: str, no_serp: bool=False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Retrieve region matrices from input cooler objects"""
    mat1 = clr1.matrix(balance=False, sparse=False).fetch(region)
    mat2 = clr2.matrix(balance=False, sparse=False).fetch(region)
    if not no_serp:
        mat2, mat1, _ = serp.serpentin_binning(
            mat2, mat1, force_symmetric=True, parallel=12, iterations=1000, threshold=10,minthreshold=2
        )
    return mat1, mat2


def parse_ucsc(region: str) -> Tuple[str, int, int]:
    """Parse a UCSC-style genomic region into a tuple of (chrom, start, end)"""
    chrom, pos = region.split(":")
    start, end = pos.split("-")
    return chrom, int(start), int(end)


def load_bed2d(path: str, region: Optional[str] = None) -> pd.DataFrame:
    """
    Load interval pairs from a bed2d file into a pandas Dataframe.
    Optionally, intervals can be filtered according to a region.
    """
    # Autodetect if there is a header
    sniffer = csv.Sniffer()
    with open(path, "r") as table:
        header = sniffer.has_header(table.read(1024))
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=list(range(6)) + [8],
        header=0 if header else None,
    )
    df.columns = [
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
        "diff_score",
    ]
    if region:
        chrom, start, end = parse_ucsc(region)
        df = df.query(
            "(chrom1 == @chrom) & (chrom2 == @chrom) & (start1 >= @start) & (start2 >= @start) & (end1 < @end) & (end2 < @end)"
        )
    return df


def get_bins(
    clr: cooler.Cooler, df: pd.DataFrame, region: Optional[str] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """Retrieve matrix bins from a dataframe of genomic intervals. Optionally,
    a region of interest can be provided to shift bins to its start.
    """
    df["ucsc1"] = df.apply(
        lambda row: f"{row.chrom1}:{row.start1}-{row.end1}", axis=1
    )
    df["ucsc2"] = df.apply(
        lambda row: f"{row.chrom2}:{row.start2}-{row.end2}", axis=1
    )
    bins1 = df["ucsc1"].apply(lambda c: clr.extent(c)[0])
    bins2 = df["ucsc2"].apply(lambda c: clr.extent(c)[0])
    if region:
        region_start_bin = clr.extent(region)[0]
        bins1 -= region_start_bin
        bins2 -= region_start_bin

    return bins1, bins2


@click.command()
@click.option(
    "--out",
    "-o",
    default=None,
    show_default=False,
    type=str,
    help="Output file where to store the heatmap plot. By default, plot interactively",
)
@click.option(
    "--matrices",
    "-m",
    help="Plot individual matrices side-by-side instead of ratio.",
    is_flag=True,
)
@click.option(
    "--no-serp",
    "-n",
    help="Use the standard log ratio instead of serpentine",
    is_flag=True,
)
@click.option(
    "--labels",
    "-l",
    help="Regions to highlight on the map. UCSC formatted strings separated by ';'.",
    type=str,
    default=None,
)
@click.option(
    "--bed2d",
    "-b",
    help="Bed2d file containing region pairs.",
    type=click.Path(exists=True),
    default=None,
)
@click.argument("cool1", type=click.Path(exists=False))
@click.argument("cool2", type=click.Path(exists=False))
@click.argument("region", type=str)
def main(
    cool1: str,
    cool2: str,
    region: str,
    out: Optional[str],
    matrices: bool = False,
    no_serp: bool = False,
    bed2d: Optional[str] = None,
    labels: Optional[str] = None,
):
    clr1 = cooler.Cooler(cool1)
    clr2 = cooler.Cooler(cool2)
    names = [basename(cool1), basename(cool2)]
    if labels:
        # Get a list of (start,end) bin tuples
        rstart = clr1.extent(region)[0]
        labels = map(clr1.extent, labels.split(";"))
        labels = [(l[0] - rstart, l[1] - rstart) for l in labels]
    if bed2d:
        # Load intervals into a dataframe
        intervals = load_bed2d(bed2d, region)
        # Compute matrix bins of each interval
        bins_x, bins_y = get_bins(clr1, intervals, region)
    if matrices:
        # Load contacts from each cooler
        m1, m2 = get_matrices(clr1, clr2, region, no_serp)
        _, ax = plt.subplots(1, 2, sharex=True, sharey=True)
        # Visualize each condition in separate panels
        for i, mat in enumerate([m1, m2]):
            ax[i].imshow(np.log10(mat), cmap="afmhot_r", rasterized=True)
            ax[i].set_title(names[i])
            # If bed2d file provided, plot each position as points
            if bed2d:
                ax[i].scatter(
                    bins_x,
                    bins_y,
                    c=intervals.diff_score,
                    cmap="PRGn",
                    edgecolors="black",
                    vmin=-0.3,
                    vmax=0.3,
                )
            if labels:
                for start, end in labels:
                    ax[i].plot(
                        [start, end], [start, end], linewidth=3, c="green"
                    )
        plt.suptitle(region)
    else:
        # Visualize a single ratio matrix
        ratio = get_ratio(clr1, clr2, region, no_serp)
        plt.imshow(
            ratio - np.nanmean(ratio),
            cmap="bwr",
            vmin=-1,
            vmax=1,
            rasterized=True,
        )
        plt.title(f"log2 {names[0]} / {names[1]}\n{region}")
        # If bed2d file provided, plot each position as points
        if bed2d:
            plt.scatter(
                bins_x,
                bins_y,
                c=intervals.diff_score,
                cmap="PRGn",
                edgecolors="black",
                vmin=-0.3,
                vmax=0.3,
            )
        plt.colorbar()
        if labels:
            for start, end in labels:
                plt.plot([start, end], [start, end], linewidth=3, c="green")
    if out is None:
        plt.show()
    else:
        plt.savefig(out)


if __name__ == "__main__":
    main()
