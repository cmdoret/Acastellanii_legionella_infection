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
    conda: '../envs/hic_processing.yaml'
    shell:
        """
        bedtools closest \
            -a <(sort -k1,1 -k2,2n {input.change}) \
            -b <(sort -k1,1 -k2,2n {input.genes}) \
        | awk 'BEGIN{{OFS="\t"}}{{print $5,$6,$7,$9,$4,$8}}' \
        > {output}
        """

# Compute percentile threshold for infection-dependent patterns
rule perc_thresh_diff_patterns:
    input: join(OUT, 'pareidolia', '{pattern}_change_infection_time.tsv')
    output: join(TMP, '{pattern}_change_thresh.txt')
    params:
        perc_thresh = 80
    shell:
        """
        tmp=$(tempfile)
        # Get number of lines (w/o header) and send sorted abs. scores to tmp file
        tot=$(
            cut -f 9 {input} \
                | tail -n+2 \
                | tr -d '-' \
                | sort -n \
                | tee $tmp \
                | wc -l
        )
        # Get row number for desired percentile
        # (n + 99) / 100 with integers is effectively ceil(n/100) with floats
        count=$(((tot * {params.perc_thresh} + 99) / 100))
        # Get score at corresponding row
        sed -n "${count}p" $tmp > {output}
        """
    

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
    conda: "../envs/viz.yaml"
    script: "../scripts/05_pattern_change_cutoff_hist.py"


# Quantify overlap of loops and borders Allowing for a margin of 1 bin
rule get_loops_overlap_borders:
    input:
        loops = join(OUT, 'pareidolia', 'loops_change_infection_time.bed'),
        borders = join(OUT, 'pareidolia', 'borders_change_infection_time.bed'),
        chroms = join(TMP, 'chrom.sizes')
    output: join(OUT, 'pareidolia', 'overlap_jaccard_loops_borders.txt')
    params:
        margin = MAX_RES
    conda: '../envs/hic_processing.yaml'
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
    script: '../scripts/05_viz_loops_overlap_borders.py'
