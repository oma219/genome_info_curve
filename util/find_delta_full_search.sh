#!/bin/bash


#Experiment 1: Try to find the cardinalities of kmer sets at varying
# k in order to be able to compute delta

#Input and Output Variables
#INPUT_FILE_ARRAY=("/home-net/home-3/oahmed6@jhu.edu/rindex_project/data/phix.fa" "/storage/oahmed/long_read_assemblies/GCA_000002125.2/genome.fna" "/storage/oahmed/general/salmonella/kuhnle/genomes/SRR1033480_contig.fa")
#INPUT_FILE_ARRAY=("/storage/oahmed/rindex_project/1OA02/salmon_concat_delta_no_canon_genome.fa")
INPUT_FILE_ARRAY=("/home-net/home-3/oahmed6@jhu.edu/rindex_project/data/phix.fa")

INPUT_FILE_NAME_ARRAY=("phix" "human" "salmonella")
#INPUT_FILE_NAME_ARRAY=("phix" "human" "salmonella")
#OUTPUT_BASE="/home-3/oahmed6@jhu.edu/rindex_project/results/exp1/"


i=-1
for INPUT_FILE in ${INPUT_FILE_ARRAY[@]}
do	
	i=$((i+1))
	INPUT_NAME="${INPUT_FILE_NAME_ARRAY[$i]}"
	OUTPUT_FILE="${OUTPUT_BASE}${INPUT_NAME}_delta_results_default_sketch.txt"
	echo $OUTPUT_FILE
	
	for kmer_size in $(seq 1 10)
	do
		curr_trial="kmer_size= ${kmer_size}"
		#echo $curr_trial >> $OUTPUT_FILE
		
		#cardcmp -k $kmer_size -S 14 --no-canon --use-cyclic-hash $INPUT_FILE >> $OUTPUT_FILE
		cardcmp -k $kmer_size -S 14 --no-canon --use-cyclic-hash $INPUT_FILE
		cardcmp -k $kmer_size -S 14 --no-canon --use-nthash  $INPUT_FILE
	done
done
