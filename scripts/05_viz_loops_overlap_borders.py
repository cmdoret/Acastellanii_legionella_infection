
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

n_loops = len(open(snakemake.input['loops'], 'r').read().split('\n'))
n_borders = len(open(snakemake.input['borders'], 'r').read().split('\n'))
stats = open(snakemake.input['jaccard'], 'r').read().split('\n')
stats = {
    k: float(v) for k, v in zip(stats[0].split('\t'), stats[1].split('\t'))
}
venn2((n_loops, n_borders, stats['n_intersections']))
plt.title(f'jaccard index: {stats["jaccard"]}')
plt.savefig(snakemake.output[0])