### This program performs Cremer-pople analysis
### Created by : Mike Limb 12/09/13  
### requires : cremer_maths_functions.py to be present as it contains the mathematical functions needed 
###	       data file with format :
### 		       #   nFr C1C2C3C4  C2C3C4C5  C3C4C5O5  C4C5O5C1  C5O5C1C2  O5C1C2C3
###    			     1  -62.576    38.474    16.080   -56.662    32.146    28.071
### usage python run_cremer_pople__puckering_analysis.py [ name of data file : eg geom_dih.txt] 
### will produce a file called Cremer_pople_puckering_coordinates.dat 
### containing .... 
### column 1 : the frame number (same as initial data file) 
### column 2 : phi angle / equatorial / azimuthal  ____ x_axis 
### column 3 : theta angle / polar 		   ____ y_axis  

from __future__ import division
import math 
import os 
#os.system("cp ~/bin/cremer_pople_tools/cremer_maths_functions.py ./")

# Define Canioincal Torsion Vectors, ie the ideal projectors for 4C1, 1,4B and OS2  

chair_proj = [60, -60, 60, -60, 60, -60]

boat_proj = [0, 60, -60, 0, 60, -60]

skew_proj = [60, -30, -30, 60, -30, -30] 

# Define Redundancy Vectors 

red4_proj= [60, 30, -30, -60, -30, 30]

red5_proj = [0, 60, 60, 0, -60, -60]

red6_proj = [60, 60, 60, 60, 60, 60]

# Calculate all the normalization factors of the ideal and redundant vectors 

import cremer_maths_functions

chair_norm = cremer_maths_functions.normalize(chair_proj)
boat_norm = cremer_maths_functions.normalize(boat_proj) 
skew_norm = cremer_maths_functions.normalize(skew_proj) 

red4_norm = cremer_maths_functions.normalize(red4_proj) 
red5_norm = cremer_maths_functions.normalize(red5_proj) 
red6_norm = cremer_maths_functions.normalize(red6_proj)

# CHECK 
# Print normalised ideal and redundant projections 
#print(chair_norm) 
#print(boat_norm) 
#print(skew_norm) 
#print(red4_norm)
#print(red5_norm) 
#print(red6_norm) 



##### Extract the Torsion vectors for the frames along the trajectories ##### 

import sys 

lines = open( sys.argv[1], "r").readlines() 

for line in range(1,len(lines)):
	torsions = lines[line].split() 
# Remove the first column (frame no) from the array so the array only contains the torsion vector 

	torsions.remove(torsions[0])
	torsions = [float(i) for i in torsions]	




##### Perform the projection coefficient calculations #####

# Calculate the chair_projection_coefficient (uncomment to check answer) 
	chair_scalar_prod = cremer_maths_functions.scalar_product(torsions, chair_proj) 
	chair_proj_coef = chair_scalar_prod / chair_norm**2
#	print(chair_proj_coef) 

# Calculate the boat_projection_coefficient (uncomment to check answer) 
	boat_scalar_prod = cremer_maths_functions.scalar_product(torsions, boat_proj)
	boat_proj_coef = boat_scalar_prod / boat_norm**2
#	print(boat_proj_coef)


# Calculate the skew(twist_boat)_projection_coefficient (uncomment to check answer)
	skew_scalar_prod = cremer_maths_functions.scalar_product(torsions, skew_proj) 
        skew_proj_coef = skew_scalar_prod / skew_norm**2
#	print(skew_proj_coef)

##### Perform the redundancy coefficient calculations #####

# Calculate the R4 coefficient (uncomment to check answer)
	red4_scalar_prod = cremer_maths_functions.scalar_product(torsions, red4_proj)
	red4_proj_coef = red4_scalar_prod / red4_norm**2
#	print(red4_proj_coef)

# Calculate the boat_projection_coefficient (uncomment to check answer) 
	red5_scalar_prod = cremer_maths_functions.scalar_product(torsions, red5_proj)
	red5_proj_coef = red5_scalar_prod / red5_norm**2
#	print(red5_proj_coef)


# Calculate the skew(twist_boat)_projection_coefficient (uncomment to check answer)
	red6_scalar_prod = cremer_maths_functions.scalar_product(torsions, red6_proj) 
        red6_proj_coef = red6_scalar_prod / red6_norm**2
#	print(red6_proj_coef)

##### Calculate Phi angle #####

	if boat_proj_coef > 0 and skew_proj_coef > 0:
		phi_radians = math.atan((skew_proj_coef * skew_norm)  / (boat_proj_coef * boat_norm) )
		phi_degrees = math.degrees(phi_radians)
		phi_degrees_abs = math.fabs(phi_degrees)
	
	elif boat_proj_coef > 0 and skew_proj_coef < 0:
		phi_radians = math.atan( (skew_proj_coef * skew_norm ) / (boat_proj_coef * boat_norm ) )
		phi_degrees = math.degrees(phi_radians)
		phi_degrees_abs = 360 - math.fabs(phi_degrees)
	
 	elif boat_proj_coef < 0 and skew_proj_coef > 0:
		phi_radians = math.atan( (skew_proj_coef * skew_norm ) / ( boat_proj_coef * boat_norm) )
		phi_degrees = math.degrees(phi_radians)
		phi_degrees_abs = 180 - math.fabs(phi_degrees)
	
	elif boat_proj_coef < 0 and skew_proj_coef < 0:
		phi_radians = math.atan( (skew_proj_coef * skew_norm) / (boat_proj_coef * boat_norm) )
		phi_degrees = math.degrees(phi_radians)
		phi_degrees_abs = 180 + math.fabs(phi_degrees)
	
	#print(phi_degrees_abs)


##### Calculate Theta angle #####

##### Calculate the possible norm_phi and then pick the smallest #####

	x0 = phi_degrees_abs - 0 
	sq_n_phi_x0 =  2*(60**2) + 2*((x0**2)) + 2*((60 - x0)**2)
	norm_phi_x0 = math.sqrt(sq_n_phi_x0)
		
	x60 = phi_degrees_abs - 60 
	sq_n_phi_x60 =  2*(60**2) + 2*((x60**2)) + 2*((60 - x60)**2)
	norm_phi_x60 = math.sqrt(sq_n_phi_x60)
	
	x120 = phi_degrees_abs - 120 
	sq_n_phi_x120 =  2*(60**2) + 2*((x120**2)) + 2*((60 - x120)**2)
	norm_phi_x120 = math.sqrt(sq_n_phi_x120)

	x180 = phi_degrees_abs - 180 
	sq_n_phi_x180=  2*(60**2) + 2*((x180**2)) + 2*((60 - x180)**2)
	norm_phi_x180 = math.sqrt(sq_n_phi_x180)

	x240 = phi_degrees_abs - 240 
	sq_n_phi_x240=  2*(60**2) + 2*((x240**2)) + 2*((60 - x240)**2)
	norm_phi_x240 = math.sqrt(sq_n_phi_x240)

	x300 = phi_degrees_abs - 300 
	sq_n_phi_x300=  2*(60**2) + 2*((x300**2)) + 2*((60 - x300)**2)
	norm_phi_x300 = math.sqrt(sq_n_phi_x300)

	norm_phi_list = [norm_phi_x0, norm_phi_x60, norm_phi_x120, norm_phi_x180, norm_phi_x240, norm_phi_x300]
	norm_phi=min(norm_phi_list) 
	#print (norm_phi)

##### Calculate the amplitude #####

	amplitude = math.sqrt(((boat_proj_coef**2)*(boat_norm**2))+((skew_proj_coef**2)*(skew_norm**2)))
	#print amplitude 

##### Do trig to get theta ######

	theta_radians = math.atan((norm_phi*chair_proj_coef) / amplitude)   
	theta_degrees = 90 + math.degrees(theta_radians) 
	
	#print(theta_degrees)

##### Put the output into a file "Cremer_pople_puckering_coordinates.dat" #####

	WFILE = open( "Cremer_pople_puckering_coordinates.dat", "a")
	print >> WFILE, "%s	 %4.9f	%4.9f" % (line, phi_degrees_abs, theta_degrees) 
	
WFILE.close()   
