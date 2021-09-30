# Rules for downloading data from the web
def lib_to_sra(wildcards):
  """
  Get SRA accession from fq path.
  """
  try:
    mask = units.fq1.str.contains(wildcards.libname).fillna(False)
    sra = units.sra[mask].values[0]
  except IndexError:
    mask = units.fq2.str.contains(wildcards.libname).fillna(False)
    sra = units.sra[mask].values[0]
  return sra

rule dl_sra:
  output: join(TMP, 'sra', '{acc}.sra')
  params:
    acc = lambda w: w.acc,
    tmp = TMP
  conda: '../envs/sra.yaml'
  singularity: 'quay.io/biocontainers/sra-tools:2.11.0--pl5262h314213e_0'
  threads: 12
  shell:
    """
    # Download SRA file
    prefetch -t {params.tmp} --max-size 100G -p -o "{output}" "{params.acc}"
    """

rule sra_to_fq:
  input: lambda w: join(TMP, 'sra', f'{lib_to_sra(w)}.sra')
  output: join('fq', '{libname}')
  message: "Getting {params.acc} into {output}"
  params:
    acc = lib_to_sra,
    tmp = lambda w: join(TMP, lib_to_sra(w))
  conda: '../envs/sra.yaml'
  singularity: 'quay.io/biocontainers/sra-tools:2.11.0--pl5262h314213e_0'
  threads: 4
  shell:
    """
    echo "TMPDIR: {params.tmp}"
    # Get library base name
    fq={output}
    trim=${{fq%_[12].fastq.gz}}
    
    # Add _1.fastq suffix if single end, otherwise, fasterq-dump adds it
    numLines=$(fastq-dump -X 1 -Z --split-spot {input} 2> /dev/null | wc -l) 
    if [ $numLines -eq 4 ]
    then
      fname="${{trim}}_1.fastq"
      echo "Extract {input} to $fname"
    else
	    fname="$trim"
      echo "Extract {input} to ${{trim}}_1.fastq and ${{trim}}_2.fastq"
    fi

    # Convert to fastq locally and compress
    mkdir -p "{params.tmp}/fq" # Workaround weird fasterq-dump deleting outdir content
    fasterq-dump -t {params.tmp} -f -e {threads} {input} -o "{params.tmp}/fq/$(basename $fname)"
    mv {params.tmp}/fq/*fastq "$(dirname {output})/"
    echo "Compress ${{trim}}*fastq"
    gzip -f ${{trim}}*fastq
    """