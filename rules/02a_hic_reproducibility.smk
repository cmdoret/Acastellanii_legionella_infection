# Assessing reproducibility of Hi-C signal between replicates of different conditions
# using HiCrep
# cmdoret, 20190930

MAXDIST = 100000

# Use two maps from replicates with best sequencing depth to select best h 
# value for future hicrep runs (PM125 vs PM54)
rule select_h_param:
    input:
        rep1 = join(OUT, 'cool', 'AT407.mcool'),
        rep2 = join(OUT, 'cool', 'AT418.mcool')
    output: join(OUT, 'hicrep', 'best_h_value.txt')
    params:
        maxdist= MAXDIST,
        res = MED_RES
    shell:
        """
        hicreppy htrain \
            {input.rep1}::/resolutions/{params.res} \
            {input.rep2}::/resolutions/{params.res} \
            -m {params.maxdist} \
            > {output}
        """


# Will be run n_lib^2 / 2 + n_lib to perform single pairwise comparisons
# Each output file will be one line: lib1 lib2 corrcoeff
# All matrices are subsampled to the same number of contacts as the
# lowest coverage sample
rule run_hicrep:
    input:
        lib1 = lambda w: join(OUT, 'cool', f'{w.library1}.mcool'),
        lib2 = lambda w: join(OUT, 'cool', f'{w.library2}.mcool'),
        h_value = join(OUT, 'hicrep', 'best_h_value.txt'),
        n_contacts = join(OUT, 'chromosight', 'target_contacts.txt')
    output: join(TMP, 'hicrep', '{library1}_{library2}_corrcoef.txt')
    params:
        max_dist = MAXDIST,
        res = MED_RES
    shell:
        """
        hicreppy scc \
                    -v $(cat {input.h_value}) \
                    -m {params.max_dist} \
                    -s $(cat {input.n_contacts}) \
                    {input.lib1}::/resolutions/{params.res} \
                    {input.lib2}::/resolutions/{params.res} \
                    > {output}
        """


# cat all corrcoeff file to get data into long format
# Load into R to generate pretty heatmap and tree
rule hicrep_matrix:
    input: expand(join(TMP, 'hicrep', '{libcombo[0]}_{libcombo[1]}_corrcoef.txt'), libcombo=it.combinations(samples.library, 2))
    output: join(OUT, 'hicrep', 'hicrep_mat.tsv')
    params:
        heatmap = join(OUT, 'plots', 'hicrep_mat.svg'),
        script = join('scripts', 'hicrep_to_mat.R'),
        samples = config['samples']
    conda: '../envs/r_env.yaml'
    shell:
        """
        # Make a single 3-column table (lib1 lib2 corrcoeff) out of all corrcoeff files
        echo -n "" > {output}
        for file in {input}; do
            tail -n+1 $file \
            | sed '/^$/d' \
            |  awk -v fname=$file 'BEGIN{{OFS="\\t"}}{{sub(".*/", "", fname);
                    split(fname,f,"_");
                    print f[1],f[2],$0}}' \
            >> {output}
        done
        # Combine all corrcoeff into a matrix
        Rscript {params.script} {output} {params.samples} {params.heatmap}
        """

