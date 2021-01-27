**********************************************
Monte Carlo Based Pathway Search Algorithm
===============================================================================

About:


Novel method used to obtain favorable trajectories using data from the exhaustive conformational search protocol while considering slow degrees of freedom (z-coordinate, inclination, azimuthal angle).  

These favorable trajectories are identified using Monte Carlo (MC) moves with a limited change in slow coordinates in each MC move. 

In our implementation: z, inclination angle, and azimuthal angle are the slow degrees of freedom.
In particular, the z-coordinates are used as a primary slow coordinate to define the pathway progression while the angles represent orthogonal, secondary slow coordinates.
The values for these variables need to be previously calculated.

This code can be run in parallel using multiple CPUs to reduce the computational time

Parameters:

		data: data from exhaustive conformational search algorithm 

		output_file: where to store the MCPS trajectories
		
		log_file: file containing errors 
		
		Indices for z, inc, azimuthal, energy in input data (inputted to MCPS initially)

			z_ind	

			inc_ind

        		az_ind

        		energy_ind

		Allowed steps for slow coordinates in each MC move
		
			z_step 
        	
			inc_step
        		
			az_step
		
		num: number of trajectories to obtain

		proc: number of processes to use to run code 
