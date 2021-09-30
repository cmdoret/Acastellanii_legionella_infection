#!/bin/env snakemake -s
from typing import List, Literal

def get_fq_units(wildcards, end: Literal[1, 2]=1) -> List[str]:
    """
    Get fastq files (units) of a particular library type of one sample 
    from the unit sheet
    """
    fqs = units.loc[
      (units.library == wildcards.library) &
      (units.libtype == wildcards.libtype),
      f"fq{end}"
    ].dropna()
    return list(fqs)

# Combine all fastq files from the same sample / library type combination
# If there is a single fastq, symlink it to spare memory
rule combine_units:
  input: lambda w: get_fq_units(w, 1)
  params:
    r2 = lambda w: get_fq_units(w, 2)
  message:
    """
    Merging :
      {input} into:  {output.r1}
      {params.r2} into:  {output.r2}
    """
  output:
    r1 = join(TMP, "reads", "{library}_{libtype}.end1.fq.gz"),
    r2 = join(TMP, "reads", "{library}_{libtype}.end2.fq.gz")
  threads: 1
  run:
    if len(input[:]) > 1:
      shell(f"cat {' '.join(input[:])} > {output['r1']}")
      if len(params['r2']):
        shell(f"cat {' '.join(params['r2'])} > {output['r2']}")
    else:
      shell(f"ln -s $PWD/{input[0]} {output['r1']}")
      if len(params['r2']):
        shell(f"ln -s $PWD/{params['r2'][0]} {output['r2']}")


# Remove small contigs from the reference
rule filter_genome:
  input: GENOME
  output: join(TMP, 'filtered_ref.fa')
  params:
    min_len = 100000
  singularity: "docker://cmdoret/seqkit:latest"
  conda: "../envs/hic_processing.yaml"
  shell: "seqkit seq -m {params.min_len} -o - {input} > {output}"
