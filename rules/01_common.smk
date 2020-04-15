#!/bin/env snakemake -s
def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    u = units.loc[
        (units.library == wildcards.library) &
        (units.libtype == wildcards.libtype)
        ,
        f'fq{wildcards.end}'
    ]
    return u.tolist()

# Combine all fastq files from the same sample / library type combination
rule combine_units:
  input: get_fastqs
  message:
    """
    Merging :
      {input} into:  {output}
    """
  output:
    join(TMP, "reads", "{library}_{libtype}.end{end}.fq.gz")
  threads: 1
  shell: "cat {input} > {output}"


# Remove small contigs from the reference
rule filter_genome:
  input: GENOME
  output: join(TMP, 'filtered_ref.fa')
  params:
    min_len = 100000
  singularity: "docker://cmdoret/seqkit:latest"
  conda: "../envs/hic_processing.yaml"
  shell: "seqkit seq -m {params.min_len} -o {output} {input}"
