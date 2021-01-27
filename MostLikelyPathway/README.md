**********************************************
Most Likely Pathway determination after performing multiple MCPS trajectories
===============================================================================

Scripts used to obtain most likely paths after multiple MCPS trajectories are identified. 

Involves 2 steps: 
1. Gridding the data and filtering trajectories passing through into groups (Filter_Trajectories.py) 
2. Evaluating a transition matrix for each group and identifying the Most Likely path from each transition matrix (Dijkstras/Find_Paths.py)

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
	
	Top and bottom z-values of grid.  Will use this grid later to select trajectories so adjust this accordingly 
        top_z
        bott_z
	
	cluster_cutoff: how to cluster trajectories.  In our case we divide the groups based on whether the azimuthal is greater/lower than 180 

	density_cutoff: we are using grids denser than this to group trajectories

	paths_file: contains all screened trajectoreis 
        paths_file_c1: cluster1 trajectories
        paths_file_c2 : cluster2 trajectories
        path_gridfle_c1 cluster1 trajectories in grid representaition ()
        path_gridfle_c2 cluster2 trajectories in grid representaition ()

	SWITCH:if you've already assigned grids to paths: SWITCH = 0, otherwise: SWITCH = 1

2. Most Likely Pathways (Find_Paths.py) 

	cluster1 transitions file
	cluster2 transitions file
	
	Indices for z, inc, azimuthal, energy in input data (inputted to MCPS initially)
	z_ind 	
        inc_in
        az_ind
        en_ind

	data: input data	

	Top and bottom z-values for algorithm.  Be careful that this is different than from grid in step 1. Will determine ending and starting z-values of Most likely path.

	step sizes used to build transition matrix.  can be different from in previous step.
	z_step
        inc_step
        az_step 

	Will output pathway files for each cluster
