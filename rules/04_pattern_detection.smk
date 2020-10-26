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
    run:
        min_contacts = np.inf
        for inp in input[:]:
            c = cooler.Cooler(inp + f"::/resolutions/{params['res']}")
            c_sum = c.info['sum']
            if c_sum < min_contacts:
                min_contacts = c_sum
            print(f"{inp}: {c_sum}")
        print(f"lowest : {min_contacts}")
        with open(output[0], 'w') as out:
            out.write(str(min_contacts))

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

rule pattern_change:
    input:
        mcools = expand( join(OUT, 'cool', '{library}.mcool'), library=samples.library),
        coords = join(OUT, 'chromosight', 'merged_contacts', '{pattern}_out.tsv')
    output: join(OUT, 'pareidolia', '{pattern}_change_infection_time.tsv')
    params:
        condition = samples.infection_time.tolist(),
        res = MAX_RES
    conda: '../envs/hic_processing.yaml'
    script: '../scripts/pattern_changes.py'

rule plot_patterns_scores:
    input: scores = expand(join(OUT, 'chromosight', '{library}', '{{pattern}}_quant.tsv'), library=samples.library)
    output: join(OUT, 'plots', '{pattern}_scores.svg')
    threads: 1
    run:
        mpl.use("Agg")
        for i, scores in enumerate(input['scores']):
            scores = pd.read_csv(scores, sep='\t')
            scores['infection_time'] = samples.infection_time.values[i]
            scores['library'] = samples.library.values[i]
            if i:
                all_scores = pd.concat([all_scores, scores])
            else:
                all_scores = scores
        sns.violinplot(data=all_scores, x='library', y='score', hue='infection_time')
        plt.ylabel(f"{wildcards['pattern']} score")
        plt.savefig(output[0])
