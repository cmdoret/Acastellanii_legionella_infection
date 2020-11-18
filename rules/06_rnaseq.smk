# We will be using using differential expression results from Li2020
# As those results are reported for the Neff genome, we first need to liftover
# Neff annotations to C3 to get a correspondance of reported genes into C3.

# Extract sheet from excel file for desired timepoint
rule get_diff_expr:
    output: join(TMP, 'diff_expr', 'de_genes.tsv')
    params:
        xls = config['rnaseq'],
        time = 8
    run: 
        time = f".{params['time']}hr"
        de = pd.read_excel(params['xls'], sheet_name=None)
        de = [v for k, v in de.items() if time in k][0]
        de = de.rename(columns={de.columns[0]: 'accession'})
        de.to_csv(output[0], sep='\t', index=False)

# parse the excel file from Li. et al 2020
# Extract the gene expression value at each time point
# Generate a new table with rows=genes and cols=time
rule get_expr_vs_time:
    output: join(OUT, 'diff_expr', 'expr_vs_time.tsv')
    params:
        xls = config['rnaseq']
    run:
        time_regex = re.compile(r'\.([0-9]+)hr')
        de = pd.read_excel(params['xls'], sheet_name=None)

        for sheet in list(de.keys()):
            # Add time column if this is a valid timepoint sheet
            try:
                timepoint = re.search(time_regex, sheet)[1]
                de[sheet]['time'] = timepoint
            # otherwise delete current sheet
            except TypeError:
                del de[sheet]
        # Concatenate all sheets into a single table
        
        de = pd.concat(de.values(), axis=0)
        (de.rename(columns={de.columns[0]: 'accession'})
            .to_csv(output[0], sep='\t', index=False)
        )

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
        # GFF files are converted to bed (using awk) and sorted 
        awk -vOFS='\t' '$3== "gene" {{
            gsub("ID=gene:","",$9);gsub(";.*", "", $9);print $1,$4,$5,$9}}' \
            {input} \
          | sort -k1,1 -k2,2n > tmp.a
        awk -vOFS='\t' '$3 == "gene" {{
            gsub("ID=","",$9);gsub(";.*", "", $9);print $1,$4,$5,$9}}' \
            {params.c3_annot} \
          | sort -k1,1 -k2,2n > tmp.b

        # The bed files are then compared to find the best overlap
        # for each gene.
        bedtools intersect -a tmp.a -b tmp.b -wao \
          | awk -vOFS='\t' '$8 != "." {{print $4,$8}}' \
          > {output}
        """

# Use the liftover to map gene identifiers from Neff v1 to C3
rule convert_identifiers:
    input:
        mapping = join(TMP, 'liftoff', 'neff_c3_gene_mapping.tsv'),
        de_genes = join(TMP, 'diff_expr', 'de_genes.tsv')
    output: join(OUT, 'diff_expr', 'de_genes.tsv')
    params:
    run:
        acc_map = pd.read_csv(input['mapping'], sep='\t', names=['neff', 'c3'])
        de_genes = pd.read_csv(input['de_genes'], sep='\t')
        de_genes = de_genes.merge(acc_map, left_on='accession', right_on='neff', how='outer')
        de_genes.drop(columns=['neff']).to_csv(output[0], sep='\t', index=False)
