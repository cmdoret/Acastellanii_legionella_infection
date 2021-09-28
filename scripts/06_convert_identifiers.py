# Use the gene mapping from liftoff to convert gene identifiers from Neff v1 to
# the new Neff assembly in the differential expression table from Li et al. 2020.
import pandas as pd


acc_map = pd.read_csv(snakemake.input['mapping'], sep='\t', names=['neff', 'c3'])
de_genes = pd.read_csv(snakemake.input['de_genes'], sep='\t')
de_genes = de_genes.merge(acc_map, left_on='accession', right_on='neff', how='outer')
(
    de_genes
    .drop(columns=['neff'])
    .to_csv(snakemake.output[0], sep='\t', index=False)
)