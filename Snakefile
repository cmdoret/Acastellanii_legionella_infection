#!/bin/env snakemake -s
# This file can be run using snakemake. It was tested on snakemake 5.3
# It orchestrates the analysis of salmonella-infected mouse macrophage.
import numpy as np
import pandas as pd
import seaborn as sns
from os.path import join
from snakemake.utils import validate
import matplotlib.pyplot as plt
import matplotlib as mpl

## LOAD CONFIG FILES
configfile: "config.yaml"
validate(config, schema='schemas/config.schema.yaml')

samples = pd.read_csv(config['samples'], sep='\t', dtype=str, comment='#').sort_values("infection_time").set_index(['library'], drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_csv(config['units'], sep='\t', dtype=str, comment='#').set_index(['library', 'unit'], drop=False)
validate(units, schema='schemas/units.schema.yaml')
# Make sure indexes are all strings
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

DATA_DIR = 'data'
IN = join(DATA_DIR, 'input')
OUT = join(DATA_DIR, 'output')
TMP = join(DATA_DIR, 'tmp')
GENOME = join(config['reference'])
MAX_RES = config['contact_maps']['max_res']
LOW_RES = config['contact_maps']['low_res']
NCPUS = config['n_cpus']

wildcard_constraints:
  library = "|".join(samples.library),
  libtype = "|".join(np.unique(units.libtype))


conda: "envs/hic_processing.yaml"

# Helper functions
include: 'scripts/mat_utils.py'
# Pipeline sub-workflows
include: 'rules/01_common.smk'
include: 'rules/02_hic_processing.smk'
include: 'rules/04_pattern_detection.smk'
include: 'rules/05_annotations_analysis.smk'

rule all:
  input:
    expand(join(OUT, 'cool', '{library}.mcool'), library=samples.library),
    expand(join(OUT, 'all_signals_{library}.bedgraph'), library=samples.library),
    join(OUT, 'plots', 'serpentine_i_u_ratio.svg'),
    expand(join(OUT, 'plots', 'coverage_hic_{library}.pdf'), library=samples.library),
    expand(join(OUT, 'plots', '{pattern}_scores.svg'), pattern=['loops', 'borders']),
    expand(join(OUT, 'plots', '{pattern}_diff_go_enrich.svg'), pattern=['loops', 'borders'])


rule aggregate_signals:
  input: 
    comp = join(OUT, 'compartments_{library}.bedgraph'),
    insul = join(OUT, 'insulation_{library}.bedgraph'),
    chroms = join(TMP, 'chrom.sizes'),
  output: join(OUT, 'all_signals_{library}.bedgraph')
  shell:
    """
    mkdir -p {TMP}/signal_tracks

    # Generate output bins as a bed file
    bedtools makewindows -w {MAX_RES} \
                         -g {input.chroms} \
                         > {OUT}/bins.bed
    
    # Intersect each signal with fixed output bins
    # Group overlapping bins and average signal for each group
    for in_file in {input.comp} {input.insul}; do
      tmp_file={TMP}/signal_tracks/$(basename $in_file)
      bedtools intersect -a $in_file -b {OUT}/bins.bed -wb |
        bedtools groupby -g 1,2,3 -c 4 -o mean 2> /dev/null |
        awk '{{print $4}}' > $tmp_file
    done
    
    echo -e "chrom\tstart\tend\tcompartment_{LOW_RES}\tlog2_insulation{LOW_RES}" > {output}
    paste {OUT}/bins.bed {TMP}/signal_tracks/*{wildcards.library}* >> {output}

    """
