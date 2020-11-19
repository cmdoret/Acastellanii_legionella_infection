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
        genes = join(IN, 'annotations', 'c3_annotations', 'Acanthamoeba_castellanii.annotations.txt')
    shell:
        """
        paste <(cut -f3-6 {params.genes}) \
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

# Compute annotation enrichment at bins from pattern changing during infection
rule go_enrich_change:
    input:
        change = join(OUT, 'pareidolia', '{pattern}_diff_genes.bed'),
        annot = join(IN, 'annotations', 'c3_annotations', 'Acanthamoeba_castellanii.annotations.txt')
    output:
        plot = join(OUT, 'plots', '{pattern}_diff_go_enrich.svg'),
        tbl = join(OUT, 'go_enrich', '{pattern}_diff_go_enrich.tsv')
    conda: "../envs/r_env.yaml"
    shell: "Rscript scripts/go_enrich.R {input.annot} {input.change} {output.plot} {output.tbl}"
