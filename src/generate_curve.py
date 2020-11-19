"""
Description: Generates curves given a list of genomes and a genomic
             repetitiveness measure such as delta, and r

Author: Omar Ahmed
"""
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

import argparse
import os
import sys
import math 
from subprocess import PIPE, Popen

def get_derivative(curr_point_stats):
	first_derv_x = (curr_point_stats[1] - curr_point_stats[0])/1
	return first_derv_x

def same_sign(x, y):
	if (x > 0 and y > 0) or (x < 0 and y < 0):
		return True
	return False

def get_optimal_ratio(start_k, start_points, mid_k,  mid_points,end_k, end_points):
	delta_value = max(max(start_points), max(mid_points), max(end_points))
	optimal_k = mid_k + mid_points.index(delta_value)

	return [optimal_k, delta_value]

def extract_result_values(stdout_split, stderr_split):
    try:
        cardinality = float(stdout_split[-1])
	elapsed_time = float(stderr_split[-3])
	system_time = float(stderr_split[-5])
	user_time = float(stderr_split[-7])
	mem_max = float(stderr_split[-1])
    except:
	print("[Error] It appears an error occurred with dashing, check log file.")
	exit(-1)

    return [cardinality, user_time, system_time, elapsed_time, mem_max]

def find_delta_binarysearch(start_k, end_k, input_genome_path, canon, log_file_name):
    
    num_calls = 0
    memo = {}
    converged = False
    mid_k = math.floor((end_k + start_k) / 2)
    iteration_num = 0
    total_user_time = 0.0
    total_system_time = 0.0
    total_elapsed_time = 0.0
    total_memory = 0.0
   
    canon_option = ""
    if not canon:
        canon_option = " --no-canon "

    while not converged:
        iteration_num += 1 
        if iteration_num > 100:
            print("Issue with Binary Search.")
	    exit(0)

        points_needed = 2
	start_points = [0 for i in range(points_needed)]
	mid_points = [0 for i in range(points_needed)]
	end_points = [0 for i in range(points_needed)]

	#Start Points
        for i in range(points_needed):
            if (start_k+i) not in memo:
                num_calls += 1
	        stream = Popen("/software/centos7/bin/time --format='user= %U system= %S elapsed= %e MemMax= %M' cardcmp -k " + str(start_k+i) + ' -S 14  --use-cyclic-hash ' +  canon_option  + input_genome_path, shell=True, stdout=PIPE, stderr=PIPE)
	        stdout, stderr = stream.communicate()
		os.popen('echo \" ' + stdout + stderr + ' \"  >> ' + log_file_name)
		cardinality, user_time, system_time, elapsed_time, mem_max  = extract_result_values(stdout.split(), stderr.split()) 
		total_user_time += user_time
		total_system_time += system_time
		total_elapsed_time += elapsed_time
		total_memory += mem_max
	
	        curr_ratio = cardinality/(start_k + i + 0.0)
                memo[(start_k+i)] = curr_ratio
            else:
                curr_ratio = memo[(start_k+i)]
	    start_points[i] = curr_ratio

        #End Points
	for i in range(points_needed):
            if (end_k+i) not in memo:
                num_calls += 1
	        stream = Popen("/software/centos7/bin/time --format='user= %U system= %S elapsed= %e MemMax= %M' cardcmp -k " + str(end_k+i) + ' -S 14  --use-cyclic-hash ' +  canon_option  +  input_genome_path, shell=True, stdout=PIPE, stderr=PIPE )
	        stdout, stderr = stream.communicate()
		os.popen('echo \" ' + stdout + stderr + ' \"  >> ' + log_file_name)
		cardinality, user_time, system_time, elapsed_time, mem_max  = extract_result_values(stdout.split(), stderr.split())
		total_user_time += user_time
		total_system_time += system_time
		total_elapsed_time += elapsed_time
		total_memory += mem_max

	        curr_ratio = cardinality/(end_k + i + 0.0)
                memo[(end_k+i)] = curr_ratio
            else:
                curr_ratio = memo[(end_k+i)]
            end_points[i] = curr_ratio

	#Mid Points
	if start_k != mid_k and mid_k != end_k:
	    for i in range(points_needed):
                if (mid_k+i) not in memo:
                    num_calls += 1
		    stream = Popen("/software/centos7/bin/time --format='user= %U system= %S elapsed= %e MemMax= %M' cardcmp -k " + str(mid_k+i) + ' -S 14  --use-cyclic-hash ' +  canon_option  +  input_genome_path, shell=True, stdout=PIPE, stderr=PIPE )
		    stdout, stderr = stream.communicate()
		    os.popen('echo \" ' + stdout + stderr + ' \"  >> ' + log_file_name)
		    cardinality, user_time, system_time, elapsed_time, mem_max  = extract_result_values(stdout.split(), stderr.split())
		    total_user_time += user_time
		    total_system_time += system_time
		    total_elapsed_time += elapsed_time
		    total_memory += mem_max

		    curr_ratio = cardinality/(mid_k + i + 0.0)
                    memo[(mid_k+i)] = curr_ratio
                else:
                    curr_ratio = memo[(mid_k+i)]
                mid_points[i] = curr_ratio

	elif start_k == mid_k:
	    mid_points = start_points
	else:
	    mid_points = end_points

	start_der = get_derivative(start_points)
	mid_der = get_derivative(mid_points)
	end_der = get_derivative(end_points)

	#Terminating Conditions
	if (start_k == mid_k) or (mid_k == end_k):
	    optimal_k, delta_value = get_optimal_ratio(start_k, start_points, mid_k,  mid_points,end_k, end_points)
	    converged = True
	else:
	    if not same_sign(start_der, mid_der):
	        end_k = mid_k
		mid_k = math.floor((end_k + start_k) / 2)
	    elif not same_sign(mid_der, end_der):
	        start_k = mid_k
		mid_k = math.floor((end_k + start_k) / 2)

    avg_memory = total_memory/num_calls
    return [optimal_k, delta_value, num_calls, total_user_time, total_system_time, total_elapsed_time, avg_memory]

def find_r(output_bwt_dir, output_genome_name, log_file_name):
    output_prefix = output_bwt_dir + "output"
    stream = Popen("/software/centos7/bin/time --format='user= %U system= %S elapsed= %e MemMax= %M' pfbwt-f64 -o " + output_prefix + " " + output_genome_name, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = stream.communicate()
    
    if log_file_name != "":
    	os.popen('echo \' ' + stderr + '\' >> ' + log_file_name)
    
    try:
	stderr_split = stderr.split()
	user_time = float(stderr_split[-7])
	system_time = float(stderr_split[-5])
	elapsed_time = float(stderr_split[-3])
	mem_max = float(stderr_split[-1])
	n = float(stderr_split[-13])
	r = float(stderr_split[-11])
    except:
 	print("[Error] It appears an error occurred with pfbwt, please check log.")
	exit(-1)
    
    return [user_time, system_time, elapsed_time, n, r, mem_max]

def write_solution(output_file, sequence, per_line=60):
    offset = 0
    while offset < len(sequence):
        nchars = min(len(sequence) - offset, per_line)
        line = sequence[offset:offset+nchars]
        offset += nchars
        output_file.write(line + "\n")

def add_reverse_comp_strand(input_file_name, output_file):

    input_file = open(input_file_name, "r")
    output_str = ""
    #Generate a complement of the string
    for line in input_file:        
        for ch in line:
            if ch == 'A':
                output_ch = 'T'
            elif ch == 'C':
                output_ch = 'G'
            elif ch == 'G':
                output_ch = 'C'
            elif ch == 'T':
                output_ch = 'A'
            elif ch == 'N':
                 output_ch = "N"
            else:
                output_ch = ""
            output_str += output_ch
        if '>' in line:
	    reverse_comp = output_str[::-1]
	    write_solution(output_file, reverse_comp)
	    output_str = ""
	    output_file.write(line + "\n")
    
    #Get the reverse complement and write it
    reverse_comp = output_str[::-1]
    write_solution(output_file, reverse_comp)
    input_file.close()

def add_forward_strand(input_file_name, output_file):

    input_file = open(input_file_name, "r")
    output_sequence = ""
    for line in input_file:
        for ch in line:
            if ch in 'ACGTN':
                output_sequence += ch
	if '>' in line:
	    write_solution(output_file, output_sequence)
	    output_sequence = ""
	    output_file.write(line + "\n")

    write_solution(output_file, output_sequence)
    input_file.close()

def preprocess_genome(input_files, file_num, output_file_name, strands_to_add, increment_size):

    output_file = open(output_file_name, "a")
    for i in range(0, increment_size):
    	#output_file.write("\n>concat_genome\n")
    	add_forward_strand(input_files[file_num+i], output_file) 
    	if strands_to_add == 1:
            #output_file.write("\n>concat_genome\n")
            add_reverse_comp_strand(input_files[file_num+i], output_file)
    output_file.close()

def generate_data(args, input_file_list):
        
    if args.log_file_name == None:
	    log_file_name = ""
    else:
	    log_file_name = args.log_file_name[0]
   
    
    for num_file in range(0, args.num_genomes[0], args.increment_size[0]):   
        
        preprocess_genome(input_file_list, num_file, args.genome_file[0], args.strands_info[0], args.increment_size[0])

        if args.measure_chosen[0] == 'r':
            user_time, system_time, elapsed_time, n, r, mem_max  = find_r(args.output_bwt_dir[0], args.genome_file[0], log_file_name)
            output_str = "\'Genomes_Added= {} user_time= {} system_time= {} elapsed_time= {} MemMax= {} n= {} r= {}\'".format(str(num_file+args.increment_size[0]), user_time, system_time, elapsed_time, mem_max, n, r)
	    os.popen('echo ' + output_str + ' >> ' + args.output_file_name[0])
    	elif args.measure_chosen[0] == 'd':
            optimal_k, delta_value, num_calls_to_find, user_time, system_time, elapsed_time, mem_max  = find_delta_binarysearch(max(1, args.kmer_seed[0]-10), args.kmer_seed[0]+10, args.genome_file[0], args.canon, log_file_name)
            args.kmer_seed = [optimal_k]
	    output_str = "\'Genomes_Added= {} user_time= {} system_time= {} elapsed_time= {} MemMax= {}  kmer_optimal= {} delta= {}\'".format(str(num_file+args.increment_size[0]), user_time, system_time, elapsed_time, mem_max, optimal_k, delta_value)
            os.popen('echo ' + output_str + ' >> ' + args.output_file_name[0])

def validate_arguments(args):

    #Series of sanity checks to make sure the options make sense
    if args.measure_chosen[0] == 'r' and args.canon == True: 
        print("[Error] --no-canon must be turned on when using r measure.")
        exit(-1)
    
    if args.measure_chosen[0] == 'r' and args.output_bwt_dir == None: 
        print("[Error] Need to specify an output directory for pfbwt when using r.")
        exit(-1)   
    
    if args.measure_chosen[0] == 'r' and args.kmer_seed != [-1]: 
        print("[Error] kmer seed only applies to calculating Delta.")
        exit(-1)
    
    if args.measure_chosen[0] == 'd' and args.kmer_seed == [-1]: 
        print("[Error] Need to specify a seed value for k to calculate Delta.")
        exit(-1)
    
    if args.measure_chosen[0] == 'd' and args.kmer_seed[0] < 1: 
        print("[Error] {} is not a valid k value to seed with.".format(args.kmer_seed[0]))
        exit(-1)

    #Make sure each file in filelist exists
    try:
        with open(args.data_list[0], "r") as input_genomes:
            input_file_list = []
            for line in input_genomes:
                line_split = line.split()
                if len(line_split) > 0 and not os.path.exists(line_split[0]):
                    print("[Error] {} cannot be found".format(line_split[0]))
                elif len(line_split) == 1 and os.path.exists(line_split[0]) and int(os.path.getsize(line_split[0])) > 0:
                    input_file_list.append(line_split[0])
    except:
        print("[Error] {} cannot be found.".format(args.data_list[0]))
        exit(-1)
    
    if len(input_file_list) == 0:
        print("[Error] No input genomes were provided.")
        exit(-1)
    
    if args.increment_size[0] < 1 or args.increment_size[0] > len(input_file_list):
        print("[Error] Not a valid increment size")
	exit(-1)

    if args.num_genomes[0] > len(input_file_list):
        print("[Error] Only {} genome(s) were provided, while {} were requested to added.".format(len(input_file_list), args.num_genomes[0]))
        exit(-1)
    
    if args.num_genomes == [-1]: #Set default to be all genomes
        args.num_genomes = [len(input_file_list)]
    
    #Try to open a genome, results, and log file if included
    try:
        output_genome_file = open(args.genome_file[0], "w")
        output_genome_file.close()
        results_file = open(args.output_file_name[0], "w")
        results_file.close()
        if args.log_file_name:
            log_file = open(args.log_file_name[0], "w")
            log_file.close()

        if args.measure_chosen[0] == 'r':
            test_output_file = open(args.output_bwt_dir[0] + "output_bwt.txt", "w")
            test_output_file.close()
            os.remove(args.output_bwt_dir[0] + "output_bwt.txt")
    except:
        print("[Error] Issue occurred when trying to open genome, results, or log file if requested.")
        exit(-1)
    
    return input_file_list
    

def arguments():
    usg = "python generate_curve.py -i data_list -m info_measure -s strands_info -g genome_file -o output_result_file -l log_output [-h] [--no-canon] [-n num_genomes] [--obwt output_bwt_dir] [--seed kmer_seed]"
    parser = argparse.ArgumentParser(description='Generate Curve', usage=usg)
    
    parser.add_argument("-i", dest= "data_list", help="Text file with name of genome to be added on each line", nargs=1)
    parser.add_argument("-m", dest="measure_chosen", help="Repetitive measure to be used, either 'r' for r or 'd' for delta", required=True, choices= ['r', 'd'], nargs= 1)
    parser.add_argument("-s", dest="strands_info", type=int, help="Strands to be included in genome, forward (0) or forward and reverse complement (1)", required= True, choices= [0, 1], nargs = 1)
    parser.add_argument("--no-canon", dest="canon", default=True, action= 'store_false', help= "When using Delta as measure, use this option to turn off kmer canonicalized [default: canonicalization]")
    parser.add_argument("-n", dest="num_genomes", type=int, help= "Number of genomes from list to use to generate curve [default: use all genomes provided]", default= [-1], nargs= 1)
    parser.add_argument("-o", dest="output_file_name", required= True, help="Name of output file", nargs= 1)
    parser.add_argument("--obwt", dest="output_bwt_dir", required= False, nargs= 1, help= "When using r measure, need to specify directory for output from pfbwt.")
    parser.add_argument("--seed", dest="kmer_seed", default= [-1], required= False, nargs= 1, type= int, help= "When using Delta measure, provide a seed value to start the delta search.")
    parser.add_argument("-l", dest="log_file_name", required=True, nargs= 1, help= "Name of log file, output from underlying used programs")
    parser.add_argument("-g", dest= "genome_file", required= True, nargs=1, help= "Name of the genome fasta file that will be created.")
    parser.add_argument("--increment", dest= "increment_size", type=int, required=True, nargs=1, help="Number of genomes that will be added before each calculation.")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = arguments()
    input_file_list = validate_arguments(args)
    generate_data(args, input_file_list)
