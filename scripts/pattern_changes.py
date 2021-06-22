
import pareidolia.hic_utils as pah

res = snakemake.params['res']
mcools = [f"{mc}::/resolutions/{res}" for mc in snakemake.input['mcools']]
bed2d = snakemake.input['coords']
conds = snakemake.params['condition']
pattern = snakemake.wildcards['pattern']

pattern_out = pah.change_detection_pipeline(
    mcools,
    conds,
    kernel=pattern,
    bed2d_file=bed2d,
    subsample=True,
    density_thresh=None,
    pearson_thresh=0.0,
    snr_thresh=0.0,
    n_cpus=8,
)

pattern_out.to_csv(snakemake.output[0], sep='\t', index=False)
