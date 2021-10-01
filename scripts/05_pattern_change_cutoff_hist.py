# Plot the histogram of pattern changes.
import matplotlib.pyplot as plt
import pandas as pd

with open(snakemake.input['thresh'], 'r') as f_thr:
	thresh = float(f_thr.read())
df = pd.read_csv(snakemake.input['change'], sep='\t')
plt.hist(df.diff_score, bins=100)
plt.axvline(thresh, c='r')
plt.axvline(-thresh, c='r')
plt.savefig(snakemake.output[0])