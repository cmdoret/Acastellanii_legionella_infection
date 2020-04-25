
#!/bin/env snakemake -s
"""
Rules from this workflow perform statistical analyses to investigate
relationships between 3D changes, annotations and gene expression.
cmdoret, 20200415
"""

# Compute enrichment of bins from detected pattern in gene annotations
rule detected_to_genes:
    input: join()
    output: join()
    script: "../scripts/"

# Compute annotation enrichment at bins from pattern changing during infection
rule diff3d_to_annots:
    input: join()
    output: join()
    script: "../scripts/"