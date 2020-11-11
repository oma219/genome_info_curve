"""
Description: Try to find delta with using the newton optimization method

Author: Omar Ahmed
"""
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import math

def get_derivatives(curr_point_stats):
	print(curr_point_stats)
	first_der_x = (curr_point_stats[1] - curr_point_stats[0])/2
	first_der_x1 = (curr_point_stats[2] - curr_point_stats[1])/2
	second_der_x = (first_der_x1 - first_der_x)/2

	return [first_der_x, second_der_x]


#---------------------------------------------------------------------
# Main Method
#---------------------------------------------------------------------

input_genome_path = "/home-net/home-3/oahmed6@jhu.edu/rindex_project/data/phix.fa"

num_iterations = 0
converged = False
curr_k = 15

while not converged:
	num_iterations += 1
	if num_iterations > 10:
		exit(0)

  	points_needed = 3
  	curr_point_stats = [0 for i in range(points_needed)]

  	for i in range(points_needed):
    		stream = os.popen('cardcmp -k ' + str(curr_k+ (2*i)) + ' -S 14  --use-cyclic-hash ' + input_genome_path + ' 2>> delta.phix.concat.run')
		prog_output = float(stream.read().split()[-1])
    		curr_ratio = prog_output/(curr_k + i + 0.0)

    		curr_point_stats[i] = curr_ratio

  	first_der, second_der = get_derivatives(curr_point_stats)
  	next_k = math.floor(curr_k - (first_der/second_der))

  	print(first_der)
  	print(second_der)
  	print(next_k)
	curr_k = next_k
