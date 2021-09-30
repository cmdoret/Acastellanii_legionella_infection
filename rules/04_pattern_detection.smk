#!/bin/env snakemake -s
"""
Rules related to loop and border detection, and how they change during infection.
cmdoret, 20200415
"""
# I merge matrices to call loops. Then compute scores of individual matrices
# in the detected coordinates
rule merge_matrices:
    input: 
        expand(
            join(OUT, 'cool', '{library}.mcool'),
            library=samples.library,
        )
    output: join(OUT, 'cool', 'all_merged.cool')
    params:
        res = MAX_RES
    conda: "../envs/hic_processing.yaml"
    threads: 8
    shell:
        """
        mcools=()
        for i in {input}; do
            echo $i
            mcools+=(${{i/%/::/resolutions/{params.res}}})
        done
        cooler merge {output} ${{mcools[@]}}
        cooler balance -p8 {output}
        """

rule detect_patterns:
    input: join(OUT, 'cool', 'all_merged.cool')
    output:
        coords = join(OUT, 'chromosight', 'merged_contacts', '{pattern}_out.tsv'),
        wins = join(OUT, 'chromosight', 'merged_contacts', '{pattern}_out.npy')
    threads: 4
    params:
        min_sep = 5 * MAX_RES,
        min_dist = lambda w: config[w.pattern]['min-dist'],
        max_dist = lambda w: config[w.pattern]['max-dist']
    conda: '../envs/hic_processing.yaml'
    shell:
        """
        out_prefix={output.coords}
        out_prefix=${{out_prefix%*.tsv}}
        chromosight detect \
            --no-plotting \
            --pattern {wildcards.pattern} \
            --win-fmt npy \
            --threads {threads} \
            --min-separation {params.min_sep} \
            --min-dist {params.min_dist} \
            --max-dist {params.max_dist} \
            --iterations 1 \
            {input} \
            $out_prefix
        """

# Matrices must be subsampled to the same coverage for loop scores to be 
# comparable. Find the lowest number of contacts among all libraries
rule find_subsampling_value:
    input:
        expand(
            join(OUT, 'cool', '{library}.mcool'),
            library=samples.library,
        )
    output: join(OUT, 'chromosight', 'target_contacts.txt')
    params:
        res = MAX_RES
    conda: '../envs/hic_processing.yaml'
    script: '../scripts/04_find_subsampling_value.py'

# Quantify loop and border scores on the individual matrices.
rule quantify_pattern_scores:
    input:
        cool = join(OUT, 'cool', '{library}.mcool'),
        coords = join(OUT, 'chromosight', 'merged_contacts', '{pattern}_out.tsv'),
        subsample = join(OUT, 'chromosight', 'target_contacts.txt')
    output:
        coords = join(OUT, 'chromosight', '{library}', '{pattern}_quant.tsv'),
        wins = join(OUT, 'chromosight', '{library}', '{pattern}_quant.npy')
    params:
        res = MAX_RES
    conda: '../envs/hic_processing.yaml'
    shell:
        """
        out_prefix={output.coords}
        out_prefix=${{out_prefix%*.tsv}}
        chromosight quantify {input.coords} \
                             {input.cool}::/resolutions/{params.res} \
                             --no-plotting \
                             --pattern {wildcards.pattern} \
                             --win-fmt npy \
                             --subsample $(cat {input.subsample}) \
                             --perc-undetected 100 \
                             --perc-zero 100 \
                             $out_prefix
        """


# Generate pseudo replicates by sampling condition-merged cools
N_REPS = 2
rule pseudo_rep_cool:
    input:
        join(OUT, 'cool', 'sub_{condition}.mcool'),
    output:
        join(OUT, 'cool', '{condition}_{i}.cool'),
    params:
        frac = 1 / N_REPS,
        res = MAX_RES
    threads: 2
    conda: '../envs/hic_processing.yaml'
    shell:
        """
        cooltools random-sample -f {params.frac} \
            {input}::/resolutions/{params.res} \
            {output}
        """


rule pattern_change:
    input:
        uni = [join(OUT, 'cool', f'uninfected_{i}.cool') for i in range(N_REPS)],
        inf = [join(OUT, 'cool', f'infected_{i}.cool') for i in range(N_REPS)],
        coords = join(OUT, 'chromosight', 'merged_contacts', '{pattern}_out.tsv')
    output: join(OUT, 'pareidolia', '{pattern}_change_infection_time.tsv')
    params:
        condition = [0 for i in range(N_REPS)] + [5 for i in range(N_REPS)]
    conda: '../envs/hic_processing.yaml'
    script: '../scripts/04_pattern_changes.py'


rule plot_patterns_scores:
    input: scores = expand(join(OUT, 'chromosight', '{library}', '{{pattern}}_quant.tsv'), library=samples.library)
    output: join(OUT, 'plots', '{pattern}_scores.svg')
    params:
        samples = samples
    threads: 1
    conda: '../envs/viz.yaml'
    script: '../scripts/04_plot_patterns_scores.py'
