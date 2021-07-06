#!/bin/env snakemake -s
# This file can be run using snakemake. It was tested on snakemake 5.3
# It orchestrates the analysis of salmonella-infected mouse macrophage.
from os.path import join
import re
import itertools as it
import numpy as np
import pandas as pd
import seaborn as sns
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
ENZ = config["enzyme"]
MAX_RES = config['contact_maps']['max_res']
MED_RES = config['contact_maps']['med_res']
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
include: 'rules/02a_hic_reproducibility.smk'
include: 'rules/04_pattern_detection.smk'
include: 'rules/05_annotations_analysis.smk'
include: 'rules/06_rnaseq.smk'

rule all:
  input:
    expand(join(OUT, 'cool', '{library}.mcool'), library=samples.library),
    join(OUT, 'plots', 'serpentine_i_u_ratio.svg'),
    #expand(join(OUT, 'plots', 'coverage_hic_{library}.pdf'), library=samples.library),
    expand(join(OUT, 'plots', '{pattern}_scores.svg'), pattern=['loops', 'borders']),
    expand(join(OUT, 'plots', '{pattern}_diff_go_enrich.svg'), pattern=['loops', 'borders']),
    join(OUT, 'hicrep', 'hicrep_mat.tsv'),
    join(OUT, 'plots', 'distance_law_infection.svg'),
    join(OUT, 'diff_expr', 'de_genes.tsv'),
    join(OUT, 'diff_expr', 'expr_vs_time.tsv'),
    expand(join(OUT, 'cool', 'sub_{condition}.mcool'), condition=['infected', 'uninfected'])


