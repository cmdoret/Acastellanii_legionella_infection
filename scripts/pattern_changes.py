
import pareidolia.hic_utils as pah
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

res = snakemake.params['res']
mcools = [f"{mc}::/resolutions/{res}" for mc in snakemake.input['mcools']]
bed2d = snakemake.input['coords']
conds = snakemake.params['condition']
graphviz = GraphvizOutput()
graphviz.output_file = 'basic.png'

with PyCallGraph(output=graphviz):
    pattern_out = pah.change_detection_pipeline(mcools, conds, bed2d_file=bed2d, subsample=True, percentile_thresh=95, n_cpus=8)

pattern_out.to_csv(snakemake.output[0], sep='\t', index=False)
