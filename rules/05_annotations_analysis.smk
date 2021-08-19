#!/bin/env snakemake -s
"""
Rules from this workflow perform statistical analyses to investigate
relationships between 3D changes, annotations and gene expression.
cmdoret, 20200415
"""

# Format infection-dependent loops into BED, with 1 anchor / line
rule loop_change_to_bed:
    input: join(OUT, 'pareidolia', 'loops_change_infection_time.tsv')
    output: temp(join(OUT, 'pareidolia', 'loops_change_infection_time.bed'))
    shell:
        """
        cat <(cut -f1-3,9 {input} | tail -n+2) \
            <(cut -f4-6,9 {input} | tail -n+2) \
        | sed 's/\t$/\t./' \
        > {output}
        """

# Format infection-dependent borders into BED, with 1 border/line
rule border_change_to_bed:
    input: join(OUT, 'pareidolia', 'borders_change_infection_time.tsv')
    output: temp(join(OUT, 'pareidolia', 'borders_change_infection_time.bed'))
    shell: "cut -f1-3,9 {input} | sed 's/\t$/\t./' | tail -n+2 > {output}"

# Format genome annotations into BED, with 1 gene/line
rule annot_to_bed:
    output: temp(join(TMP, 'genes.bed'))
    params:
        genes = join(IN, 'annotations', 'c3_annotations', 'Acanthamoeba_castellanii_C3.annotations.txt')
    shell:
        """
        paste <(cut -f4-6 {params.genes}) \
              <(cut -f1 {params.genes}) \
        | sed 's/\t$/\t./' \
        | tail -n+2 \
        > {output}
        """

# Report the closest gene to each infection-dependent pattern
rule bed_closest_genes:
    input:
        change = join(OUT, 'pareidolia', '{pattern}_change_infection_time.bed'),
        genes = join(TMP, 'genes.bed')
    output: join(OUT, 'pareidolia', '{pattern}_diff_genes.bed')
    shell:
        """
        bedtools closest \
            -a <(sort -k1,1 -k2,2n {input.change}) \
            -b <(sort -k1,1 -k2,2n {input.genes}) \
        | awk 'BEGIN{{OFS="\t"}}{{print $5,$6,$7,$9,$4,$8}}' \
        > {output}
        """

# Define a threshold for infection-dependent patterns
rule perc_thresh_diff_patterns:
    input: join(OUT, 'pareidolia', '{pattern}_change_infection_time.tsv')
    output: join(TMP, '{pattern}_change_thresh.txt')
    params:
        perc_thresh = 80
    run:
        df = pd.read_csv(input[0], sep='\t')
        thresh = np.percentile(df.diff_score.abs(), params['perc_thresh'])
        with open(output[0], 'w') as f_out:
            f_out.write(str(thresh))
    

# Compute annotation enrichment at bins from pattern changing during infection
rule go_enrich_change:
    input:
        change = join(OUT, 'pareidolia', '{pattern}_diff_genes.bed'),
        annot = join(IN, 'annotations', 'c3_annotations', 'Acanthamoeba_castellanii_C3.annotations.txt'),
        thresh = join(TMP, '{pattern}_change_thresh.txt')
    output:
        plot = join(OUT, 'plots', '{pattern}_diff_go_enrich.svg'),
        tbl = join(OUT, 'go_enrich', '{pattern}_diff_go_enrich.tsv')
    conda: "../envs/r_env.yaml"
    shell: "Rscript scripts/go_enrich.R {input.annot} {input.change} {input.thresh} {output.plot} {output.tbl}"

# Show histogram of pattern change along with cutoff
rule pattern_change_cutoff_hist:
    input:
        change = join(OUT, 'pareidolia', '{pattern}_change_infection_time.tsv'),
        thresh = join(TMP, '{pattern}_change_thresh.txt')
    output: join(OUT, 'plots', '{pattern}_diff_cutoff_hist.svg')
    run:
        with open(input['thresh'], 'r') as f_thr:
            thresh = float(f_thr.read())
        df = pd.read_csv(input[0], sep='\t')
        plt.hist(df.diff_score, bins=100)
        plt.axvline(thresh, c='r')
        plt.axvline(-thresh, c='r')
        plt.savefig(output[0])


# Quantify overlap of loops and borders Allowing for a margin of 1 bin
rule get_loops_overlap_borders:
    input:
        loops = join(OUT, 'pareidolia', 'loops_change_infection_time.bed'),
        borders = join(OUT, 'pareidolia', 'borders_change_infection_time.bed'),
        chroms = join(TMP, 'chrom.sizes')
    output: join(OUT, 'pareidolia', 'overlap_jaccard_loops_borders.txt')
    params:
        margin = MAX_RES
    shell:
        """
        bedtools jaccard \
            -a <(sort -k1,1 -k2,2n {input.borders} \
                | bedtools slop -i - -g {input.chroms} -b {params.margin}) \
            -b <(sort -k1,1 -k2,2n {input.loops} \
                | bedtools slop -i - -g {input.chroms} -b {params.margin}) \
            > {output}
        """


rule viz_loops_overlap_borders:
    input:
        jaccard = join(OUT, 'pareidolia', 'overlap_jaccard_loops_borders.txt'),
        loops = join(OUT, 'pareidolia', 'loops_change_infection_time.bed'),
        borders = join(OUT, 'pareidolia', 'borders_change_infection_time.bed')
    output: join(OUT, 'plots', 'venn_loops_borders_overlap.svg')
    run:
        n_loops = len(open(input['loops'], 'r').read().split('\n'))
        n_borders = len(open(input['borders'], 'r').read().split('\n'))
        stats = open(input['jaccard'], 'r').read().split('\n')
        stats = {k: float(v) for k, v in zip(stats[0].split('\t'), stats[1].split('\t'))}
        from matplotlib_venn import venn2
        venn2((n_loops, n_borders, stats['n_intersections']))
        plt.title(f'jaccard index: {stats["jaccard"]}')
        plt.savefig(output[0])
