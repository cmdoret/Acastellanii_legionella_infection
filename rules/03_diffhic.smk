#!/bin/env snakemake -s

DTMP = join(TMP, 'diffhic')

# Combine single end bam files into paired end file and fix mate info.
rule prep_diffhic_bam:
  input:
    end1 = join(TMP, 'bam', '{library}_hic.end1.bam'),
    end2 = join(TMP, 'bam', '{library}_hic.end2.bam')
  output: temporary(join(DTMP, 'bam','{library}.bam'))
  threads: 1
  shell: "python scripts/stitch_bams.py -1 {input.end1} -2 {input.end2} -o {output}"

# Prepare input pairs format for diffHic
rule prep_diffhic_pairs:
  input: expand(join(DTMP, 'bam', '{library}.bam'), library=samples.library)
  output: expand(join(DTMP, 'pairs', '{library}.h5'), library=samples.library)
  params:
    bam_dir = join(DTMP, 'bam'),
    pairs_dir = join(DTMP, 'pairs')
  shell: "Rscript scripts/prep_diffhic_input.R -i {params.bam_dir} -o {params.pairs_dir}"

# Compute differential domain boundaries between infected and uninfected samples
rule diff_domain_boundaries:
  input: expand(join(DTMP, 'pairs', '{library}.h5'), library=samples.library)
  output: join(OUT, 'diffhic', 'diff_domain_boundaries.txt')
  params:
    pairs_dir = join(DTMP, 'pairs'),
  shell: "Rscript scripts/diffhic_contacts.R -i {params.pairs_dir} -o {output}"

# Get significant regions (pvalue < 0.01) into BED format
rule diffhic_to_bed:
  input: join(OUT, 'diffhic', 'diff_domain_boundaries.txt')
  output: join(OUT, 'diffhic', 'sig_diff_domain_boundaries.bed')
  params:
    threshold = 0.01
  run:
    dom = pd.read_csv(input[0], sep='\t')
    sig = dom.loc[dom.PValue< params['threshold'], ['seqnames', 'start', 'end', 'logFC']] 
    sig.to_csv(output[0], sep='\t', index=False)

rule serpentine_binning:
  input:
    a = join(OUT, 'cool', "AT337.mcool"),
    b = join(OUT, 'cool', "PM106.mcool")
  output: join(OUT, 'plots', 'serpentine_i_u_ratio.svg')
  params:
    serp_res = LOW_RES
  shell:
    """
    python scripts/serpentine_analysis.py \
      {input.a}::/resolutions/{params.serp_res} \
      {input.b}::/resolutions/{params.serp_res} \
      {output}
    """
