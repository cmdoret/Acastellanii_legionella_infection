#!/bin/env snakemake -s

def make_units_dict():
  """
  Used to work around the incompatibility between google cloud and functions as input.
  Builds a dictionary of paths for librairies and libtypes and prepends google storage 
  bucket name to each path.
  """
  # path_dict = {libtype: {library: {fq1: [path1, ...], fq2: [path1, ...]}}}
  path_dict = {}
  for t in np.unique(units.libtype):
    path_dict[t] = {}
    for s in samples.library:
      path_dict[t][s] = {}
      for e in ["fq1", "fq2"]:
        fqs = units.loc[(units.library == s) & (units.libtype == t), e].dropna().values.tolist()
        path_dict[t][s][e] = fqs
  return path_dict

units_dict = make_units_dict()

# Combine all fastq files from the same sample / library type combination
rule combine_units:
  input: lambda w: units_dict[f'{w.libtype}'][f'{w.library}'][f'fq{w.end}']
  message:
    """
    Merging :
      {input} into:  {output}
    """
  output:
    join(TMP, "reads", "{library}_{libtype}.end{end}.fq.gz")
  threads: 12
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
