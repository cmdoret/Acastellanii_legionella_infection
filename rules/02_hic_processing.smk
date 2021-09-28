#!/bin/env snakemake -s

rule bt2_index:
  input: join(TMP, 'filtered_ref.fa')
  output: touch(join(TMP, 'genome.bt2.done'))
  params:
    idx = join(TMP, 'genome')
  singularity: "docker://koszullab/hicstuff:3.1.0"
  conda: "../envs/hic_processing.yaml"
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
  singularity: "docker://cmdoret/hicstuff:3.1.0"
  conda: "../envs/hic_processing.yaml"
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
  conda: "../envs/hic_processing.yaml"
  shell: "samtools merge -n -O BAM -@ {threads} {output} {input}"


## 00 Generate Hi-C pairs files
rule generate_pairs:
  input:
    bam1 = join(TMP, 'bam', "{library}_hic.end1.bam"),
    bam2 = join(TMP, 'bam', "{library}_hic.end2.bam")
  output:
    hicdir = directory(join(TMP, 'hicstuff', '{library}')),
    pairs = join(TMP, 'pairs', '{library}.pairs')
  params:
    enz = ENZ,
    idx = join(TMP, 'genome')
  threads: 1
  singularity: "docker://koszullab/hicstuff:3.1.0"
  conda: "../envs/hic_processing.yaml"
  log: "logs/hicstuff/{library}.log"
  shell:
    """
    hicstuff pipeline --force \
                      -e {params.enz} \
                      -g {params.idx} \
                      -o {output.hicdir} \
                      -S bam \
                      --filter \
                      -P {wildcards.library} \
                      -nD \
                      {input.bam1} \
                      {input.bam2} 2> {log}

    cp {output.hicdir}/tmp/{wildcards.library}.valid_idx_pcrfree.pairs {output.pairs}
    """


# Convert pairs to cool
rule pairs_to_cool:
  input:
    chroms = join(TMP, "chrom.sizes"),
    pairs = join(TMP, 'pairs', '{library}.pairs')
  output: join(OUT, 'cool', '{library}.cool')
  params:
    res = MAX_RES
  conda: "../envs/hic_processing.yaml"
  shell:
    """
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
                       {input.chroms}:{params.res} \
                       {input.pairs} {output}
                       
    """


rule digest_genome:
  input: join(TMP, 'filtered_ref.fa')
  output:
    tigs = join(TMP, 'digest', 'info_contigs.txt'),
    frags = join(TMP, 'digest', 'fragments_list.txt')
  params:
    enz = ENZ
  conda: "../envs/hic_processing.yaml"
  shell:
    """
    hicstuff digest --force \
			        -e {params.enz} \
              -o $(dirname {output.frags}) \
              {input}
    """


# 00: Generate chrom sizes file
rule chrom_sizes:
  input:
    ref = join(TMP, 'filtered_ref.fa'),
    tigs = join(TMP, 'digest', 'info_contigs.txt')
  output: join(TMP, "chrom.sizes")
  shell:
    """
    tail -n +2 {input.tigs} \
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
  conda: "../envs/hic_processing.yaml"
  shell:
    """
		cooler zoomify -r {params.max_res},{params.med_res},{params.low_res} \
				       --balance \
					   --balance-args '--mad-max 10' \
					   -n {threads} \
					   {input.cool}
    """

# Merge ends of Hi-C bam files and sort by coord
rule merge_sort_bam_ends:
  input: expand(join(TMP, 'bam', '{{library}}_hic.{end}.bam'), end=['end1', 'end2'])
  output: temporary(join(TMP, 'bam', '{library}_hic.merged.bam'))
  threads: NCPUS
  conda: "../envs/hic_processing.yaml"
  shell:
    """
	samtools merge -n -@ {threads} -O BAM - {input} \
	  | samtools sort -@ {threads} -O BAM -o {output}
    samtools index -@ {threads} {output}
    """


# 04: Visualise Hi-C coverage along genome for each library
rule plot_hic_coverage:
  input: join(TMP, 'bam', '{library}_hic.end1.bam')
  output:
      plot = join(OUT, 'plots', 'coverage_hic_{library}.pdf'),
	  text = join(OUT, 'cov_hic', 'coverage_hic_{library}.bedgraph')
  params:
    win_size = 100000,
    win_stride = 10000
  threads: 12
  conda: "../envs/hic_processing.yaml"
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
  conda: '../envs/hic_processing.yaml'
  shell:
    """
    python scripts/02_serpentine_analysis.py \
      {input.a}::/resolutions/{params.serp_res} \
      {input.b}::/resolutions/{params.serp_res} \
      {output}
    """

rule compute_genomic_distance_law:
  input: 
    pairs = join(TMP, 'pairs', '{library}.pairs'),
    frags = join(TMP, 'digest', 'fragments_list.txt')
  output:
    tbl = join(TMP, 'distance_law', '{library}_ps.tsv'),
    plt = temp(join(TMP, 'distance_law', '{library}_ps.svg'))
  conda: '../envs/hic_processing.yaml'
  shell:
    """
    hicstuff distancelaw -a \
                         -i 1000 \
                         -f {input.frags} \
                         --pairs {input.pairs} \
                         -O {output.tbl} \
                         -o {output.plt}
    """


rule plot_distance_law:
  input: expand(join(TMP, 'distance_law', '{library}_ps.tsv'), library=samples.library)
  output: join(OUT, 'plots', 'distance_law_infection.svg')
  params:
    names = ','.join(samples.library.values.tolist()),
    inputs = ','.join([join(TMP, 'distance_law', f'{lib}_ps.tsv') for lib in samples.library])
  conda: '../envs/hic_processing.yaml'
  shell: "hicstuff distancelaw -a -l {params.names} -o {output} --dist-tbl='{params.inputs}'"


# Compute total contacts by condition for subsampling
rule get_merged_contacts:
  input:
    lambda w: expand(
      join(OUT, 'cool', '{library}.cool'),
      library=samples.loc[samples.condition==w.condition, 'library']
    )
  output: join(TMP, '{condition}_contacts.txt')
  conda: '../envs/hic_processing.yaml'
  shell:
    """
    declare -i contacts
    for cl in {input}
    do
      curr=$(cooler info $cl \
        | grep sum \
        | sed 's/[^0-9]*\([0-9]\+\)$/\1/')
      contacts=$((contacts+curr))
    done
    contacts > {output}
    """


# Retrieve the smallest number of contacts between conditions for subsampling
rule get_lowest_contacts:
  input: expand(join(TMP, '{condition}_contacts.txt'), condition=['infected', 'uninfected'])
  output: join(TMP, 'lowest_condition_contacts.txt')
  shell: "cat {input} | sort -k1,1n | head -1 > {output}"


# Generate condition-merged cool files, subsampled at the same coverage for both conditions.
rule condition_merged_cools:
  input:
    cool = lambda w: expand(
      join(OUT, 'cool', '{library}.cool'),
      library=samples.loc[samples.condition==w.condition, 'library']
    ),
    target_contacts = join(TMP, 'lowest_condition_contacts.txt')
  output: join(OUT, 'cool', 'sub_{condition}.mcool')
  threads: 3
  params:
      out_cool = lambda w: join(OUT, 'cool', f'sub_{w.condition}.cool'),
      max_res = MAX_RES,
      med_res = MED_RES,
      low_res = LOW_RES
  conda: "../envs/hic_processing.yaml"
  shell:
    """
    COUNT=$(head -n1 {input.target_contacts})
    cooler merge {params.out_cool}_tmp {input.cool}
    # Skips subsampling for the condition with lowest coverage
    cooltools random-sample --count "$COUNT" {params.out_cool}_tmp {params.out_cool} \
      || cp {params.out_cool}_tmp {params.out_cool}
		cooler zoomify -r {params.max_res},{params.med_res},{params.low_res} \
		  --balance \
		  --balance-args '--mad-max 10' \
		  -n {threads} \
		  {params.out_cool}
    rm {params.out_cool} {params.out_cool}_tmp
    """
