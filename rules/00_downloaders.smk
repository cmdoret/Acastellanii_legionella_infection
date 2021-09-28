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


rule sra_dl_fq:
  message: "Getting {params.acc} into {output}"
  output: join('fq', '{libname}')
  params:
    acc = lib_to_sra
  conda: '../envs/sra.yaml'
  singularity: 'quay.io/biocontainers/sra-tools:2.11.0--pl5262h314213e_0'
  threads: 12
  shell:
    """
    # Check if output files are already present in case rule was called by other mate
    if [ ! -f {output} ]
    then
      # Download SRA file
      prefetch -t ./fq --max-size 100G -p -o "./fq/{params.acc}.sra" "{params.acc}"
      
      # Get library base name
      fq={output}
      trim=${{fq%_[12].fastq.gz}}
      
      # Add _1.fastq suffix if single end, otherwise, fasterq-dump adds it
      numLines=$(fastq-dump -X 1 -Z --split-spot ./fq/{params.acc}.sra | wc -l) 
      if [ $numLines -eq 4 ]
      then
        trim="${{trim}}_1.fastq"
        echo "SRA download to ${{trim}}"
      else
        echo "SRA download to ${{trim}}_1.fastq and ${{trim}}_2.fastq"
      fi

      # Convert to fastq locally and compress
      fasterq-dump -t ./fq -f -e {threads} "./fq/{params.acc}.sra" -o $trim
      rm -f "./fq/{params.acc}.sra"
      gzip -f ${{trim}}*fastq
    fi
    """