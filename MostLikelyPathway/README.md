**********************************************
Most Likely Pathway Identification 
===============================================================================

Scripts used to obtain most likely paths after multiple MCPS trajectories are identified. 

Involves 2 steps: 

1. Gridding the trajectory data and filtering trajectories passing through high density grids into groups according to a cutoff (Filter_Trajectories.py) 
2. Evaluating a transition matrix for each group and applying Dijkstra's algorithm to each transition matrix to identify the most Likely path  (Dijkstras/Find_Paths.py)

Parameters for each step:

1. Filtering trajectories (Filter_Trajectories.py)

		data_file: exhaustive search data

		transition_file: all transitions from MCPS trajectories

		num_transitions: number of trajectories to use (use number after which convergence is reached)

		Indices for z, inc, azimuthal, energy in input data (inputted to MCPS initially)
	
			z_index	
	
			inc_index
	
	        	az_index
	
	        	en_index
	
		length of the grid in z, inclination, azimuthal dimensions
	
			z_step 
		
	        	inc_step
	
	        	az_step
		
		Top and bottom z-values of grid.  Will use this grid to select trajectories so adjust this accordingly 
	
		        top_z
		        bott_z
		
		density_cutoff: we are using grids denser than this to group trajectories.
	
		cluster_cutoff: how to group trajectories.  In our case we divide the groups based on whether the azimuthal is greater/lower than a cutoff
	
		The points on each trajectory need to be represented as grids. SWITCH controls whether you need to assign these grids or you have already done so.
		
			SWITCH:if you've already assigned grids: SWITCH = 0, otherwise: SWITCH = 1
	
		paths_file: all screened trajectories file
	        paths_file_c1: cluster1 trajectories file
	        paths_file_c2 : cluster2 trajectories file
	
	     	Final output: filtered trajectories for each group 

2. Most Likely Pathways (Find_Paths.py) 

		cluster1 transitions file

		cluster2 transitions file

		Indices for z, inc, azimuthal, energy in input data 

			z_ind 	
        		
			inc_in
        		
			az_ind

        		en_ind

		data: input data	

		Top and bottom z-values to apply Dijkstra's algorithm.  Be careful that this is different than from grid in step 1. Will determine ending and starting z-values of Most likely path.

			top_z

			bott_z

		Step sizes used to build transition matrix. 
			
			z_step
        		
			inc_step
        		
			az_step 

		Final output:  most likely pathway for each cluster
