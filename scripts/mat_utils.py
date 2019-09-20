# Utilities to compute stats from sparse matrices in bedgraph2 format
# 20190430, cmdoret

import pandas as pd
import numpy as np
from hicstuff import io as hio
from hicstuff import hicstuff as hcs
import cooler
import h5py
import bioframe
from cooltools import eigdecomp


def coverage_to_bedgraph(mat_file, cov_file, binsize=5000):
    """
    Compute coverage along genome from an input matrix in bedgraph2 format
    and writes the results to a bedgraph file.

    Parameters
    ----------
    mat_file : str
        Path to the input matrix, in bedgraph2D format.
    cov_file : str
        Path to the ouput bedgraph file to generate.
    binsize : int
        Bin size in which to compute coverage, in basepairs.
    
    Examples
    --------
    coverage_to_bedgraph("mat.bg2", "coverage.bg", 1000)
    """

    coo_mat, frags = hio.load_bedgraph2d(mat_file)
    cov_sum = coo_mat.sum(axis=1)


def compartment_to_bedgraph(cool, bedgraph):
    """
    Computes eigen vectors from a cooler file, correlate and ranks them according
    to gene density per bin in hg19 and writes a bedgraph file.

    Parameters
    ----------
    cool : str
        Path to the input cooler file. Can also be the URI in an mcool file. For
        example: ex.mcool::/resolutions/640000
    bedgraph : str
        Path to the output bedgraph file with compartment info.
    """
    # Retrieve cooler file and extract bin table
    c = cooler.Cooler(cool)
    bins = c.bins()[:]
    # Fetch and compute gene coverage per bin and make a new bins table with
    # gene_count column
    genecov = bioframe.tools.frac_gene_coverage(bins, "hg19")
    # Compute 3 eigen vectors and rank + correlate with gene density
    cis_vals, cis_eigs = eigdecomp.cooler_cis_eig(
        c, bins=genecov, n_eigs=3, phasing_track_col="gene_count"
    )

    bg = cis_eigs.loc[:, ["chrom", "start", "end", "E1"]]
    bg.to_csv(bedgraph, sep="\t", header=None, index=False, na_rep="nan")


def add_chromosome_arms(cool, centro):
    """
    Writes a 'chrom_arm' column to the "bins' table of a cool file.

    Parameters
    ----------
    cool : str
        Path to a cool file.
    centro : str
        Path to a BED file with centromeres regions.
    """
    c = cooler.Cooler(cool)
    # Read centromere file and assign each in to a centromere
    centro = pd.read_csv(centro, sep='\t', header=None, names=['chrom', 'centro_start', 'centro_end'])
    centro['centro_mid'] = (centro['centro_start'] + centro['centro_end']) // 2
    arm = c.bins()[:].merge(centro, how='left', suffixes=['b', 'c']).loc[:, ['start', 'centro_mid']]
    # Define the chromosome arm of each bin (B if bin start > middle of centromere, else A)
    arm['arm'] = np.where(arm['start'] > arm['centro_mid'], 'B', 'A')
    # Contigs that have no centromeres are given NaN as chromosome arm
    arm.loc[np.isnan(arm.centro_mid), 'arm'] = np.nan
    # Add new feature to cooler
    hio.add_cool_column(c, arm['arm'], 'chrom_arm', 'bins', dtype=h5py.special_dtype(vlen=str))
