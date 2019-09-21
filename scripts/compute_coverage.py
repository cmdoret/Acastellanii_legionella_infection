# Compute and plot coverage in rolling windows along the genome of input bam file.
# cmdoret, 20190920

import os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pysam as ps
import click


def parse(sam, chromlist, res):
    """
    Parse input sam or bam file and yield chromosomes one by one along with a
    rolling window mean of coverage.

    Parameters
    ----------
    sam : pysam.Samfile
        A pysam file handle to an alignment file.
    chromlist : list of str
        A list of chromosome names to include in the analysis.

    Returns
    -------
    generator of str, int, pandas.Series
        For each element in the generator, there are 3 values: The chromosome
        name, its length and an array of rolling coverage values.
    """
    for chromo, length in zip(sam.references, sam.lengths):
        if chromo in chromlist:
            depths = np.zeros(length + 1)
            for base in sam.pileup(chromo, setpper="all"):
                try:
                    depths[base.reference_pos] += base.nsegments
                except AttributeError:
                    depths[base.pos] += len(base.pileups)
            df = pd.DataFrame(depths, index=np.arange(length + 1), columns=["depth"])
            yield chromo, length, df.rolling(window=res, center=True).mean()


@click.command()
@click.option(
    "--res",
    "-r",
    default=10000,
    help="Size of windows in which to compute coverage, in basepairs.",
)
@click.option(
    "--skip", "-s", default=1000, help="Stride between windows, in basepairs."
)
@click.option(
    "--name",
    "-n",
    default="",
    help="Name of the sample (plot title). Base name of input file by default",
)
@click.option(
    "--blacklist",
    "-b",
    default="",
    help="Exclude those chromosomes from the plot. List of comma-separated chromosome names.",
)
@click.option(
    "--whitelist",
    "-w",
    default="",
    help="Only include those chromosomes in the plot. List of comma-separated chromosome names.",
)
@click.option(
    "--out",
    "-o",
    default="",
    help="Output file where to write the plot. If not provided, the plot is shown interactively",
    type=click.Path(exists=False),
)
@click.argument("bam", type=click.Path(exists=True))
def covplot(bam, out, res, skip, name, blacklist, whitelist):
    click.echo("Visualise read coverage in rolling windows from a bam file.")
    sns.set_style("white")
    scale = 1000
    sam = ps.Samfile(bam)
    blacklist = blacklist.split(",")
    whitelist = whitelist.split(",")
    chromlist = []
    if len(whitelist[0]):
        chromlist = whitelist
    else:
        chromlist = list(sam.references)
        if len(blacklist[0]):
            for chrom in blacklist:
                chromlist.remove(chrom)

    with sns.color_palette("husl", sam.nreferences):
        min_count, max_count = 0, 0
        offset, chrom_id = np.zeros(len(chromlist)+1), 0
        for chrom, length, counts in parse(sam, chromlist, res):
            plt.plot(
                (counts.index.values[::skip] + offset[chrom_id]) / scale,
                counts[counts.columns[0]].values[::skip],
                ".",
            )
            highest = np.max(counts.iloc[::skip, 0])
            lowest = np.min(counts.iloc[::skip, 0])
            if lowest < min_count:
                min_count = lowest
            if highest > max_count:
                max_count = highest

            plt.axvline(offset[chrom_id] / scale)
            offset[chrom_id+1] = offset[chrom_id] + length
            chrom_id += 1
    for n, chrom in enumerate(chromlist):
            plt.text(((offset[n+1] - offset[n]) / 2 + offset[n]) / scale, 1.05 * max_count , chrom)
    plt.xlabel("kb")
    # plt.legend()
    plt.gca().set_ylim([min_count, 1.1 * max_count])
    plt.ylabel(f"coverage ({res}-bp averaged)")
    plt.gca().set_xlim([0, offset[-1] / scale])
    if len(name) == 0:
        name = os.path.splitext(os.path.basename(bam))[0]
    plt.title(name)
    if len(out):
        plt.savefig(out)
    else:
        plt.show()


if __name__ == "__main__":
    covplot()
