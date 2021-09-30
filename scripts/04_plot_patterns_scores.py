# Visualize the distribubtion of pattern scores for both conditions as violin
# plots
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

samples = snakemake.params['samples']
mpl.use("Agg")
for i, scores in enumerate(snakemake.input['scores']):
    scores = pd.read_csv(scores, sep='\t')
    scores['infection_time'] = samples.infection_time.values[i]
    scores['library'] = samples.library.values[i]
    if i:
        all_scores = pd.concat([all_scores, scores])
    else:
        all_scores = scores
sns.violinplot(data=all_scores, x='library', y='score', hue='infection_time')
plt.ylabel(f"{snakemake.wildcards['pattern']} score")
plt.savefig(snakemake.output[0])