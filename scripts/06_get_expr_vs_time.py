# Parse the table from Li et al. 2020 to retrieve expression versus time for all genes.
import pandas as pd
import re


time_regex = re.compile(r'\.([0-9]+)hr')
de = pd.read_excel(snakemake.params['xls'], sheet_name=None)

for sheet in list(de.keys()):
    # Add time column if this is a valid timepoint sheet
    try:
        timepoint = re.search(time_regex, sheet)[1]
        de[sheet]['time'] = timepoint
    # otherwise delete current sheet
    except TypeError:
        del de[sheet]
# Concatenate all sheets into a single table

de = pd.concat(de.values(), axis=0)
(de
    .rename(columns={de.columns[0]: 'accession'})
    .to_csv(snakemake.output[0], sep='\t', index=False)
)