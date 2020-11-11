"""
Description:

Author: Omar Ahmed

"""

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import math

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


#---------------------------------------------------------------------
# Main Method
#---------------------------------------------------------------------

input_genome_path = "/home-net/home-3/oahmed6@jhu.edu/rindex_project/data/phix.fa"

converged = False
start_k = 1
end_k = 100
mid_k = math.floor((end_k + start_k) / 2)
iteration_num = 0

while not converged:
	iteration_num += 1
	print("Iteration Num = " +str(iteration_num))
	if iteration_num > 10:
		exit(0)

	points_needed = 2
	start_points = [0 for i in range(points_needed)]
	mid_points = [0 for i in range(points_needed)]
	end_points = [0 for i in range(points_needed)]
	
	#Start Points
        for i in range(points_needed):
		stream = os.popen('cardcmp -k ' + str(start_k+i) + ' -S 14  --use-cyclic-hash ' + input_genome_path + ' 2>> delta.phix.concat.run')
		prog_output = float(stream.read().split()[-1])
		curr_ratio = prog_output/(start_k + i + 0.0)
		start_points[i] = curr_ratio

        #End Points
	for i in range(points_needed):
		stream = os.popen('cardcmp -k ' + str(end_k+i) + ' -S 14  --use-cyclic-hash ' + input_genome_path + ' 2>> delta.phix.concat.run')
		prog_output = float(stream.read().split()[-1])
		curr_ratio = prog_output/(end_k + i + 0.0)
		end_points[i] = curr_ratio

	#Mid Points 
	if start_k != mid_k and mid_k != end_k:
		for i in range(points_needed):
			stream = os.popen('cardcmp -k ' + str(mid_k+i) + ' -S 14  --use-cyclic-hash ' + input_genome_path + ' 2>> delta.phix.concat.run')
			prog_output = float(stream.read().split()[-1]) 
			curr_ratio = prog_output/(mid_k + i + 0.0)
			mid_points[i] = curr_ratio
	elif start_k == mid_k:
		mid_points = start_points
	else:
		mid_points = end_points
	

	start_der = get_derivative(start_points)
	mid_der = get_derivative(mid_points)
	end_der = get_derivative(end_points)
 	
	print("Starting Points:")
	print(start_k)
	print(mid_k)
	print(end_k)
	
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

	print("Derivatives: ")
	print(start_der)
	print(mid_der)
	print(end_der)

print(optimal_k)
print(delta_value)
