####################################################
# Name: exp_1.smk
# Description: Creates an ordering of genomes, and 
#              iteratively add in genomes, and 
#              compute the number of runs and phrases
#              Lempel-Ziev parse.
#
# Date: 10/31/22
####################################################


####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_all_genomes_in_dataset_exp1(wildcards):
    """ Returns a list of all the genomes """
    input_files = []
    for data_file in os.listdir(f"data/dataset_1/"):
        if data_file.endswith(".fna"):
            input_files.append(f"data/dataset_1/{data_file}")
    return input_files

####################################################
# Section 3: Rules needed for this experiment type
####################################################

# Section 3.1: Build the list of genomes, the filelist will be
# overall ordering used for every reference. The second rule
# builds a refernce using the first n genomes in the filelist.

rule build_genome_list_exp1:
    input:
        get_all_genomes_in_dataset_exp1
    output:
        "exp1_filelist/genome_filelist.txt"
    run:
        with open(output[0], "w") as out_fd:
            for path in input:
                out_fd.write(f"{path}\n")

rule build_reference_file_with_x_genomes_exp1:
    input:
        "exp1_filelist/genome_filelist.txt"
    output:
        "exp1_ref_files/ref_with_{num}_genomes.fna"
    shell:
        """
        genomes_included=0
        while IFS= read -r path; do
            cat $path >> {output}

            genomes_included=$((genomes_included+1))
            if (( $genomes_included == {wildcards.num})); then
                break
            fi
        done < {input}
        """

# Section 3.2: Builds the BWT and reports the number of runs
# and the time it took for a specific reference.

rule build_bwt_for_x_genomes_exp1:
    input:
        "exp1_ref_files/ref_with_{num}_genomes.fna"
    output:
        "exp1_workspace/{num}_genomes_r/ref_with_{num}_genomes.fna.ri",
        "exp1_workspace/{num}_genomes_r/ref_with_{num}_genomes.fna.ri.log"
    shell:
        """
        cp {input} exp1_workspace/{wildcards.num}_genomes_r/ref_with_{wildcards.num}_genomes.fna

        cd {r_dir} 
        ri-buildfasta {base_dir}/exp1_workspace/{wildcards.num}_genomes_r/ref_with_{wildcards.num}_genomes.fna > {base_dir}/{output[1]}
        cd {base_dir}
        """

rule extract_r_data_for_single_reference_exp1:
    input:
        "exp1_workspace/{num}_genomes_r/ref_with_{num}_genomes.fna.ri.log"
    output:
        "exp1_raw_data/r_data/ref_{num}_genomes.txt"
    shell:
        """
        time=$(grep "construction time:" {input} | awk '{{print $NF}}')
        r=$(grep "equal-letter" {input}| awk '{{print $NF}}')
        printf "%d %.3f\n" "$r" "$time" > {output}
        """

# Section 3.3: Builds the LZ77 parse for a specific reference, and
# reports the number of phrases in that parse.

rule build_lz77_for_x_genomes_exp1:
    input:
        "exp1_ref_files/ref_with_{num}_genomes.fna"
    output:
        "exp1_workspace/{num}_genomes_z/ref_with_{num}_genomes.fna.parse",
        "exp1_workspace/{num}_genomes_z/ref_with_{num}_genomes.fna.parse.log.1",
        "exp1_workspace/{num}_genomes_z/ref_with_{num}_genomes.fna.parse.log.2"
    shell:
        """
        cp {input} exp1_workspace/{wildcards.num}_genomes_z/ref_with_{wildcards.num}_genomes.fna
        inputfasta="{base_dir}/exp1_workspace/{wildcards.num}_genomes_z/ref_with_{wildcards.num}_genomes.fna"

        cd {lz77_dir}
        _deps/bigbwt-build/newscanNT.x $inputfasta -w 10 -p 100 -f > {base_dir}/{output[1]}
        src/lz_77_test64 $inputfasta > {base_dir}/{output[2]}
        cd {base_dir}
        """

rule extract_z_data_for_single_reference_exp1:
    input:
        "exp1_workspace/{num}_genomes_z/ref_with_{num}_genomes.fna.parse.log.1",
        "exp1_workspace/{num}_genomes_z/ref_with_{num}_genomes.fna.parse.log.2"
    output:
        "exp1_raw_data/z_data/ref_{num}_genomes.txt"
    shell:
        """
        echo "hello" > {output}
        z=$(grep "phrase number" {input[1]} | awk '{{print $NF}}')
        time1=$(grep "Elapsed time" {input[1]} | tail -n1 | awk '{{print $NF}}')
        time2=$(grep "Elapsed time:"  {input[0]} | awk '{{print $4}}')
        totaltime=$(echo "$time1 + $time2" | bc)

        printf "%d %.3f\n" "$z" "$totaltime" > {output}
        """

# Section 3.4: Combines all the r- and z-data into 
# a single csv file.

rule gather_all_raw_data_for_csv_exp1:
    input:
        expand("exp1_raw_data/z_data/ref_{num}_genomes.txt", num=genome_breaks),
        expand("exp1_raw_data/r_data/ref_{num}_genomes.txt", num=genome_breaks)
    output:
        "exp1_final_data/output.csv"
    shell:
        """
        printf "numgenomes,r,rtime,z,ztime\n" > {output}
        for num in {genome_breaks}; do
            r_file="exp1_raw_data/r_data/ref_${{num}}_genomes.txt"
            z_file="exp1_raw_data/z_data/ref_${{num}}_genomes.txt"

            r=$(head -n1 $r_file | awk '{{print $1}}')
            z=$(head -n1 $z_file | awk '{{print $1}}')

            rtime=$(head -n1 $r_file | awk '{{print $2}}')
            ztime=$(head -n1 $z_file | awk '{{print $2}}')
            
            printf "%d,%d,%.3f,%d,%.3f\n" $num $r $rtime $z $ztime >> {output}
        done
        """






# ./_deps/bigbwt-build/newscanNT.x ../data/yeast.fasta -w 10 -p 100 -f
# ./src/lz_77_test64 ../data/yeast.fasta



# "exp1_ref_files/ref_with_{num}_genomes.fna"

#     # Expected Output: numgenomes,r,z,r-time,z-time
# "exp1_ref_files/ref_with_{num}_genomes.fna"

#     output:
#         "exp7_combined_refs/dataset_{num}/combined_ref_all.fna.ri"
#     shell:
#         """
#         cd {r_dir} 
#         ri-buildfasta {base_dir}/{input}
#         cd {base_dir}
        # """

# [INFO] 19:50:03 - Message: Window size set to:  10
# [INFO] 19:50:03 - Message: Computing LZ77
# [INFO] 19:50:03 - Message: Computing SA, LCP, and DA of dictionary
# [INFO] 19:50:23 - Message: Elapsed time (s):  19.919
# [INFO] 19:50:23 - Message: Computing ISA of dictionary
# [INFO] 19:50:26 - Message: Elapsed time (s):  2.49284
# [INFO] 19:50:26 - Message: Computing RMQ over LCP of dictionary
# [INFO] 19:50:27 - Message: Elapsed time (s):  1.10907
# [INFO] 19:50:27 - Message: Computing co-lex DA of dictionary
# [INFO] 19:50:36 - Message: Elapsed time (s):  9.2818
# [INFO] 19:50:36 - Message: Computing SA of the parsing
# [INFO] 19:50:36 - Message: Elapsed time (s):  0.106763
# [INFO] 19:50:36 - Message: Computing ISA of the parsing
# [INFO] 19:50:36 - Message: Elapsed time (s):  0.00305001
# [INFO] 19:50:36 - Message: Computing LCP of the parsing
# [INFO] 19:50:36 - Message: Elapsed time (s):  0.0037841
# [INFO] 19:50:36 - Message: Computing RMQ over LCP of the parsing
# [INFO] 19:50:36 - Message: Elapsed time (s):  0.00703258
# [INFO] 19:50:36 - Message: Computing b_p
# [INFO] 19:50:36 - Message: Elapsed time (s):  0.0260999
# [INFO] 19:50:36 - Message: Computing b_bwt, b_pps, and M of the parsing
# [INFO] 19:50:56 - Message: Elapsed time (s):  19.5079
# [INFO] 19:50:56 - Message: Computing S_LCP_T
# [INFO] 19:50:56 - Message: Elapsed time (s):  0.579659
# [INFO] 19:50:56 - Message: create lz77
# [INFO] 19:50:57 - Message: tree size:  489594
# [INFO] 19:50:57 - Message: after all the support data structures, RAM:  3752152688
# [INFO] 19:52:46 - Message: LZ77 construction:  2996075236
# [INFO] 19:52:46 - Message: phrase number:  1472523
# [INFO] 19:52:46 - Message: LZ77 construction time (s):  109.393
# [INFO] 19:52:46 - Message: Elapsed time (s):  109.605
# [INFO] 19:52:46 - Message: LZ77 construction complete
# [INFO] 19:52:46 - Message: Elapsed time (s):  163.566
# malloc_count ### exiting, total: 14,236,478,432, peak: 2,996,075,685, current: 973,770