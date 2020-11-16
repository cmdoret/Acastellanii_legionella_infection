# We will be using using differential expression results from Li2020
# As those results are reported for the Neff genome, we first need to liftover
# Neff annotations to C3 to get a correspondance of reported genes into C3.

# Note: the current assembly url seems broken in genomepy 
# we will use it when fixed, for now, the genome path must be provided in config
rule download_neff:
    output: join(TMP, 'liftover', 'neff_v1.fa')
    params:
        accession = "Acastellanii.strNEFF_v1"
    conda: '../envs/genomepy.yaml'
    shell:
        """
        genomepy install {params.accession}
        """


# Liftover annotations from the published Neff (v1) genome
# to the current C3 assembly
rule liftover_annotations:
    output: join(TMP, 'liftoff', 'neff_c3_liftover.gff')
    params:
        neff_fa = config['neff']['genome'],
        neff_annot = config['neff']['annot'],
        c3_fa = GENOME
    conda: '../envs/liftoff.yaml'
    threads: NCPUS
    shell:
        """
        liftoff -g {params.neff_annot} \
            -p {threads} \
            -o {output} \
            {params.c3_fa} \
            {params.neff_fa}
        """

# intersect liftover coordinates with de novo C3 annotations
# to get mapping from Neff to C3 identifiers
rule map_neff_c3_identifiers:
    input: join(TMP, 'liftoff', 'neff_c3_liftover.gff')
    output: join(TMP, 'liftoff', 'neff_c3_gene_mapping.tsv')
    params:
        c3_annot = config['annot']
    shell:
        """
        # GFF files are converted to bed (using awk) and sent
        # to bedtools using process substitution
        # The bed file is then sorted and only the best overlap
        # for each gene is retained.
        bedtools intersect \
            -a <( awk '$3== "gene" {{gsub("ID=gene:","",$9);print $1$4,$5,$9}}' {input} ) \
            -b <( awk '$3 == "gene" {{gsub("ID=","",$9);print $1,$4,$5,$9}}' {params.c3_annot} ) \
            -wo \
          | sort -k4,8 -k9,9rn \
          | cut -f4,8 \
          | uniq > {output}
        """

# Use the liftover to map gene identifiers from Neff to C3
rule convert_identifiers:
    input: join(TMP, 'liftoff', 'neff_c3_gene_mapping.tsv')
    output: join(OUT, 'diff_expr', 'de_genes.tsv')
    params:
        de_genes = config['de_genes']
    run:
        acc_map = pd.read_csv(input[0], sep='\t')
        de_genes = pd.read_csv(params['de_genes'], sep='\t')
        de_genes.merge(acc_map, right_index=True, left_index=True, how='left')
        de_genes.to_csv(output[0], sep='\t')