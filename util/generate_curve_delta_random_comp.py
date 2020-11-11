#Author: Omar Ahmed
#Description: Reads in genome files and incrementally concantenates them and calculates
#some statistic over the genome file that can be used to generate a curve

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import math
import random
import copy
from copy import deepcopy

def add_forward_strand_same_base(input_file_name, base_comps, output_same_base_file):
    
    base_dict = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    forward_str = ""
    while sum(base_comps) > 0:
        base_index = random.choice([0, 1, 2, 3])

        while base_comps[base_index] == 0:
            base_index = random.choice([0, 1, 2, 3])

        base_comps[base_index] -= 1
        curr_base = base_dict[base_index]
        forward_str += curr_base
    
    output_same_base_file.write(forward_str)

def add_forward_strand_random_base(input_file_name, base_comps,  output_random_base_file):
    
    base_dict = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    forward_str = ""
    genome_chars_left = sum(base_comps)
    
    while genome_chars_left > 0:

        base_index = random.choice([0, 1, 2, 3])
        curr_base = base_dict[base_index]
        forward_str += curr_base

        genome_chars_left -= 1

    output_random_base_file.write(forward_str)
    

def add_forward_strand(input_file_name, output_file):
    
    base_comps = [0, 0, 0, 0]
    input_file = open(input_file_name, "r")
    for line in input_file:
        output_line = ""
        for ch in line:
            if ch in 'ACGTN':
                output_line += ch
            
            if ch == 'A':
                base_comps[0] += 1
            elif ch == 'C':
                base_comps[1] += 1
            elif ch == 'G':
                base_comps[2] += 1
            elif ch == 'T':
                base_comps[3] += 1
        
        output_file.write(output_line)
    input_file.close()
    return base_comps

def preprocess_genome(input_files, file_num, output_file_name, output_same_base, output_random_base):

    output_file = open(output_file_name, "a")
    output_file.write("\n>salmon_concat_genome\n")
    base_comps = add_forward_strand(input_files[file_num], output_file)
    output_file.close()

    output_same_base_file = open(output_same_base, "a")
    output_same_base_file.write("\n>salmon_concat_genome\n")
    add_forward_strand_same_base(input_files[file_num], deepcopy(base_comps), output_same_base_file)
    output_same_base_file.close()
    
    output_random_base_file = open(output_random_base, "a")
    output_random_base_file.write("\n>salmon_concat_genome\n")
    add_forward_strand_random_base(input_files[file_num], deepcopy(base_comps), output_random_base_file)
    output_random_base_file.close()

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

def find_delta_binarysearch(start_k, end_k, input_genome_path):
    
    num_calls = 0
    memo = {}
    converged = False
    mid_k = math.floor((end_k + start_k) / 2)
    iteration_num = 0

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
	        stream = os.popen('cardcmp -k ' + str(start_k+i) + ' -S 14  --use-cyclic-hash ' + input_genome_path + ' 2>> delta.phix.concat.run')
	        prog_output = float(stream.read().split()[-1]) 
	        curr_ratio = prog_output/(start_k + i + 0.0)
                memo[(start_k+i)] = curr_ratio
            else:
                curr_ratio = memo[(start_k+i)]
	    start_points[i] = curr_ratio

        #End Points
	for i in range(points_needed):
            if (end_k+i) not in memo:
                num_calls += 1
	        stream = os.popen('cardcmp -k ' + str(end_k+i) + ' -S 14  --use-cyclic-hash ' + input_genome_path + ' 2>> delta.phix.concat.run')
	        prog_output = float(stream.read().split()[-1])
	        curr_ratio = prog_output/(end_k + i + 0.0)
                memo[(end_k+i)] = curr_ratio
            else:
                curr_ratio = memo[(end_k+i)]
            end_points[i] = curr_ratio

	#Mid Points
	if start_k != mid_k and mid_k != end_k:
	    for i in range(points_needed):
                if (mid_k+i) not in memo:
                    num_calls += 1
		    stream = os.popen('cardcmp -k ' + str(mid_k+i) + ' -S 14  --use-cyclic-hash ' + input_genome_path + ' 2>> delta.phix.concat.run')
		    prog_output = float(stream.read().split()[-1])
		    curr_ratio = prog_output/(mid_k + i + 0.0)
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


    return [optimal_k, delta_value, num_calls]

#-----------------------------------------------------------------------
# Main Method
#-----------------------------------------------------------------------


#Parameters for the program
num_of_files = 100
kmer_mid = 14
kmer_mid_2 = 14
kmer_mid_3 = 14
output_genome_name = "salmon_concat_genome_delta_random_exp.fa"
output_same_base = "salmon_concat_genome_delta_samebasecomp.fa"
output_random_base = "salmon_concat_genome_delta_randombasecomp.fa"
output_stats_name = "salmon_concat_delta_results_random.txt"

#Remove concatenated genome file from previous run, just so we don't add to it by mistake
if os.path.exists(output_genome_name):
    os.remove(output_genome_name)
if os.path.exists(output_same_base):
    os.remove(output_same_base)
if os.path.exists(output_random_base):
    os.remove(output_random_base)

#Grab a list of genome files
genome_dir = '\'/storage/oahmed/general/salmonella/kuhnle/genomes/\''
stream = os.popen('find ' + genome_dir + ' -name \'*.fa\' ')
input_file_list = stream.read().split()

#Open file to write down delta statistics
output_stats_file = open(output_stats_name, "w")
output_stats_file.write("num of genome files = " + str(num_of_files) + "\n")
output_stats_file.close()

#Remove any empty genomes
for file_name in input_file_list:
    size = int(os.path.getsize(file_name))
    if size == 0:
        input_file_list.remove(file_name)

for i in range(num_of_files):

    #Go through next genome and add it ...
    preprocess_genome(input_file_list, i, output_genome_name, output_same_base, output_random_base)

    #Run binary search to get optimal k and delta
    optimal_k, delta_value, num_calls_to_find = find_delta_binarysearch(max(1, kmer_mid-10), kmer_mid+10, output_genome_name)
    kmer_mid = optimal_k

    optimal_k_2, delta_value_2, num_calls_to_find_2 = find_delta_binarysearch(max(1, kmer_mid_2-10), kmer_mid_2+10, output_same_base)
    kmer_mid_2 = optimal_k_2

    optimal_k_3, delta_value_3, num_calls_to_find_3 = find_delta_binarysearch(max(1, kmer_mid_3-10), kmer_mid_3+10, output_random_base)
    kmer_mid_3 = optimal_k_3

    output_str = "\'Salmonella  --> Iteration_Num= " + str(i) + ' num_calls= '  + str(num_calls_to_find)  + ' kmer_size= '+ str(optimal_k) + ' Delta= ' + str(delta_value)
    output_str += "\nSame_Base   --> Iteration_Num= " + str(i) + ' num_calls= '  + str(num_calls_to_find_2)  + ' kmer_size= '+ str(optimal_k_2) + ' Delta= ' + str(delta_value_2)
    output_str += "\nRandom_Base --> Iteration_Num= " + str(i) + ' num_calls= '  + str(num_calls_to_find_3)  + ' kmer_size= '+ str(optimal_k_3) + ' Delta= ' + str(delta_value_3) + "\n\'"
    os.popen('echo ' + output_str + ' >> ' + output_stats_name)

#output_stats_file.close()
