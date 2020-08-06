#!/bin/env snakemake -s


rule bt2_index:
  input: join(TMP, 'filtered_ref.fa')
  output: touch(join(TMP, 'genome.bt2.done'))
  params:
    idx = join(TMP, 'genome')
  singularity: "docker://koszullab/hicstuff:latest"
  shell: "bowtie2-build {input} {params.idx}"

# Make splits from Hi-C fastq files to speed up mapping. [between 1 and 999]
N_SPLITS = 4
split_names = [f'part_{s:03}' for s in range(1, N_SPLITS + 1)] #base names of split files

rule split_hic_fastq:
  input: join(TMP, 'reads', '{library}_hic.{end}.fq.gz')
  output: expand(join(TMP, 'split_reads', '{{library}}_hic_{{end}}', '{{library}}_hic.{{end}}.{split}.fq.gz'), split=split_names)
  params:
    n_splits = N_SPLITS,
    split_dir = lambda w: join(TMP, 'split_reads', f"{w.library}_hic_{w.end}")
  message: "Splitting {wildcards.library}_{wildcards.end} into {params.n_splits} split fastq"
  singularity: "docker://cmdoret/seqkit:latest"
  conda: "../envs/hic_processing.yaml"
  threads: 1
  shell:
     """
     mkdir -p {params.split_dir}
     # 100 split fastqs will be created with name pattern 00000.fq - 000100.fq
     seqkit split2 -p {params.n_splits} \
                   -w 0 \
                   -f \
                   -1 {input} \
                   -O {params.split_dir}
     """

# Alignment of a single fastq split from a Hi-C library
rule split_align_hic:
  input:
    index_flag = join(TMP, 'genome.bt2.done'),
    fq = join(TMP, 'split_reads', '{library}_hic_{end}', '{library}_hic.{end}.{split}.fq.gz'),
  output: join(TMP, 'split_reads', '{library}_hic_{end}', '{library}_hic.{end}.{split}.bam')
  params:
    index = join(TMP, 'genome'),
    bt2_presets = config['params']['bowtie2']
  threads: 12
  singularity: "docker://cmdoret/hicstuff:latest"
  shell:
    """
    bowtie2 {params.bt2_presets} \
            -p {threads} \
            -x {params.index} \
            -U {input.fq} | 
      samtools sort -@ {threads} -n -O BAM -o {output}
    """

# Merge splits from individual mapping jobs into one bam file per library end
rule merge_split_alignments:
  input:
    expand(
      join(TMP, 'split_reads', '{{library}}_hic_{{end}}', '{{library}}_hic.{{end}}.{split}.bam'),
      split=split_names
    )
  output: join(TMP, 'bam', '{library}_hic.{end}.bam')
  threads: 12
  singularity: "docker://biocontainers/samtools:v1.7.0_cv4"
  shell: "samtools merge -n -O BAM -@ {threads} {output} {input}"

## 00 Generate Hi-C pairs files
rule generate_pairs_cool:
  input:
    bam1 = join(TMP, 'bam', "{library}_hic.end1.bam"),
    bam2 = join(TMP, 'bam', "{library}_hic.end2.bam")
  output:
    hicdir = directory(join(TMP, 'hicstuff', '{library}')),
    pairs = join(TMP, 'pairs', '{library}.pairs'),
    cool = join(OUT, 'cool', '{library}.cool')
  params:
    res = config['contact_maps']['max_res'],
    idx = join(TMP, 'genome')
  threads: 1
  singularity: "docker://koszullab/hicstuff:latest"
  shell:
    """
    hicstuff pipeline --force \
				      -e {params.res} \
                      -g {params.idx} \
                      -o {output.hicdir} \
                      -S bam \
                      -P {wildcards.library} \
                      -nD \
                      -M cool \
                      {input.bam1} \
                      {input.bam2}

    cp {output.hicdir}/tmp/{wildcards.library}.valid_idx_pcrfree.pairs {output.pairs}
    cp {output.hicdir}/{wildcards.library}.cool {output.cool}
    """

# 00: Generate chrom sizes file
rule chrom_sizes:
  input: join(TMP, 'filtered_ref.fa')
  output: join(TMP, "chrom.sizes")
  params:
    genome = join(TMP, 'filtered_ref.fa'),
    res = MAX_RES,
    digest_dir = join(TMP, 'digest')
  singularity: "docker://cmdoret/hicstuff:latest"
  shell:
    """
    hicstuff digest --force \
			        -e {params.res} \
                    -o {params.digest_dir} \
                    {input}

    tail -n +2 {params.digest_dir}/info_contigs.txt \
      | awk -vOFS='\t' '{{print $1,$2}}' \
      > {output}
    """

# 02: Generate multiple resolutions and normalize mcool files 
# (increments of factor 2)
rule zoomify_normalize_cool:
  input:
    cool = join(OUT, 'cool', '{library}.cool')
  output: join(OUT, 'cool', '{library}.mcool')
  threads: 3
  params:
      max_res = MAX_RES,
      med_res = MED_RES,
      low_res = LOW_RES
  singularity: "docker://cmdoret/cooler:0.8.5"
  shell:
    """
		cooler zoomify -r {params.max_res},{params.med_res},{params.low_res} \
				       --balance \
					   --balance-args '--mad-max 10' \
					   -n {threads} \
					   {input.cool}
    """


# 03c: Compute insulation scores along the matrix (TAD boundaries) the diamond 
# window is set to 50x50 pixels
rule insulation_score:
  input: join(OUT, 'cool', '{library}.mcool')
  output: join(OUT, 'insulation_{library}.bedgraph')
  params:
    win_size_bp = 10 * LOW_RES,
    res = LOW_RES
  singularity: "docker://cmdoret/cooler:0.8.5"
  shell: 
    """
    cooltools diamond-insulation {input}::/resolutions/{params.res} \
                                        {params.win_size_bp} |
      tail -n +2 |
      awk -vOFS="\t" '{{print $1,$2,$3,$5}}' > {output}
    """

# Merge ends of Hi-C bam files and sort by coord
rule merge_sort_bam_ends:
  input: expand(join(TMP, 'bam', '{{library}}_hic.{end}.bam'), end=['end1', 'end2'])
  output: temporary(join(TMP, 'bam', '{library}_hic.merged.bam'))
  threads: NCPUS
  shell:
    """
	samtools merge -n -@ {threads} -O BAM - {input} \
	  | samtools sort -@ {threads} -O BAM -o {output}
    samtools index -@ {threads} {output}
    """
# 04: Visualise Hi-C coverage along genome for each library
rule plot_hic_coverage:
  input: join(TMP, 'bam', '{library}_hic.merged.bam')
  output:
      plot = join(OUT, 'plots', 'coverage_hic_{library}.pdf'),
	  text = join(OUT, 'cov_hic', 'coverage_hic_{library}.bedgraph')
  params:
    win_size = 10000,
    win_stride = 100
  shell:
    """
	tinycov covplot \
      -N \
      -n {wildcards.library} \
      -s {params.win_stride} \
      -r {params.win_size} \
	  -t {output.text} \
      -o {output.plot} \
      {input}
    """

rule serpentine_binning:
  input:
    a = join(OUT, 'cool', "AT421.mcool"),
    b = join(OUT, 'cool', "AT420.mcool")
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
