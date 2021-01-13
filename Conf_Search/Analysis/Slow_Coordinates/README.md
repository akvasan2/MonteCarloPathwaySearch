**********************************************
Evaluate Pair Interaction Energy (PIE) for each minimized pose
===============================================================================
About:
	
	Scripts used to analyze the slow coordinates of our minimized pose data.  For antibiotic permeation case, we use the z-axis, inclination angle, and azimuthal angle

	Very important to obtain this data to use in the MCPS algorithm

Required files:
	
	Structure (i.e. psf) file stored as ../../Input_Files/system.psf
	
	All minimized pose dcd files stored in ../../Minimization/minimize_output/dcdsstart_*.dcd 

Important Parameters to specify in slow_coor.tcl: 
	
	head_name and tail_name: head and tail of reference vector for calculation.
	output_file: output file name
	num_dcds: number of minimized pose data dcds

How to run:
	
	vmd slow_coor.tcl 

Output:
	
	slow_coor.dat: columns are z-coor, inclination, azimuthal data for each pose 
