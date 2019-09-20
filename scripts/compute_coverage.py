# Compute and plot coverage in rolling windows along the genome of input bam file.
# cmdoret, 20190920

import sys
import os.path
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pysam as ps
import click

def parse(sam, chromlist):
    """
    Parse input sam or bam file and yield chromosomes one by one along with a rolling window
    mean of coverage.

    Parameters
    ----------
    sam : pysam.Samfile
        A pysam file handle to an alignment file.
    chromlist : list of str
        A list of chromosome names to include in the analysis.

    Returns
    -------
    generator of str, int, pandas.Series
        For each element in the generator, there are 3 values: The chromosome name, its length
        and an array of rolling coverage values.
    """
    for chromo, length in zip(sam.references, sam.lengths):
        depths = np.zeros(length + 1)
        for base in sam.pileup(chromo, setpper='all'):
            try:
                depths[base.reference_pos] += base.nsegments
            except AttributeError:
                depths[base.pos] += len(base.pileups)
        df = pd.DataFrame(depths, index=np.arange(length + 1), columns=['depth'])
        yield chromo, length, df.rolling(window=res).mean()

@click.command()
@click.option('--res', default=10000, help='Size of windows in which to compute coverage, in basepairs.')
@click.option('--skip', default=1000, help='Stride between windows, in basepairs.')
@click.option('--name', default='', help='Name of the sample (plot title). Base name of input file by default',
@click.option('--blacklist', default='', help='Exclude those chromosomes from the plot. List of comma-separated chromosome names.',
@click.option('--whitelist', default='', help='Only include those chromosomes in the plot. List of comma-separated chromosome names.')
@click.argument('bam')
@click.argument('out')
def covplot(bam, out, res, skip, name, blacklist, whitelist):
    click.echo("Visualise read coverage in rolling windows in a bam file.")
    sns.set_style('white')
    sam = ps.Samfile(bam)
    blacklist = blacklist.split(',')
    whitelist = whitelist.split(',')
    chromlist = []
    if len(whitelist):
        chromlist = whitelist
    else:
        chromlist = list(sam.refrences)
        for chrom in blacklist:
            chromlist.remove(chrom)

    with sns.color_palette('husl', sam.nreferences):
        offset = 0
        for chrom, length, counts in parse(sam, chromlist):
            plt.plot(counts.index.values[::skip] + offset, counts[counts.columns[0]].values[::skip], '.') 
            offset += length

            medians[chromo] = np.median(counts[counts.columns[0]])
            all_depths.extend(counts[counts.columns[0]][::skip])
    
    plt.xlab('bp')
    plt.legend()
    plt.ylab(f'coverage ({res}-bp averaged)')
    plt.gca().set_xlim([0, offset+res+1])
    if len(name) == 0:
        name = os.path.splitext(os.path.basename(bam))[0]
    plt.title(name)
    plt.savefig(out)

if __name__ == "main__":
    covplot()
