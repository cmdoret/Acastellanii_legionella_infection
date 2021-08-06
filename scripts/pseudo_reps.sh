# Commands used to make pseudo-replicates.
# cmdoret, 20210805

# Merge fastqs for each condition separately
for end in 1 2; do
	cat fq/{AT407,AT438,PM106,PM106_reseq,AT418,AT420}.end${end}.fq.gz \
		> fq/infected_merged.end${end}.fq.gz
	cat fq/{AT337,AT408,AT419,AT421}.end${end}.fq.gz \
		> fq/uninfected_merged.end${end}.fq.gz
done

# Random sample of merged fastq with different seeds for each pseudo-replicate
mkdir -p fq_pseudo
for end in 1 2; do
	for rep in {1..3}; do
		seqkit sample -p 0.3 -s $rep fq/uninfected_merged.end${end}.fq.gz \
			-j6 -o fq_pseudo/pu${rep}.end${end}.fq.gz
		seqkit sample -p 0.3 -s $rep fq/infected_merged.end${end}.fq.gz \
			-j6 -o fq_pseudo/pi${rep}.end${end}.fq.gz
	done
done
