# Parse data from Li et al. 2020 to extract differential expression values at
# desired timepoint from the corresponding excel sheet.
import pandas as pd

time = f".{snakemake.params['time']}hr"
de = pd.read_excel(snakemake.params['xls'], sheet_name=None)
de = [v for k, v in de.items() if time in k][0]
de = de.rename(columns={de.columns[0]: 'accession'})
de.to_csv(snakemake.output[0], sep='\t', index=False)