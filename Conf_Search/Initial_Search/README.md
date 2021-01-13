**********************************************
Exhaustive search for antibiotic poses in conformational space 
===============================================================================
About:
	
	Scripts used to exhaustively search over the translational and rotational space for an antibiotic.
	
	This entails running 2 protocols:
		
		1. Generation of multiple drug orientations by aligning a vector of the drug to fibonacci sphere and doing self rotations of the drug
		2. Placing the drug in all possible positions within the protein at points along a grid
			Removes any clases between drug and protein
	
	A third protocol is called to remove ring piercings

Required files:
	
	Structure and coordiante files (i.e. psf, pdb) files stored as ../Input_Files/system.psf and ../Input_Files/system.pdb

Important Parameters to specify in config.tcl: 
	
	1. Generation of multiple orientations:
	
		num: number of fibonacci points used to generate angle values
		
		head_name and tail_name: head and tail of reference vector for rotation of antibiotic.

		self_rot: angle specified for self rotation of the angle 

	2. Translation of drug to points on a grid:

		pocket_residues: residues of the protein which are used to choose a center for a grid

		del_x_pos, del_x_neg, del_y_pos, ..., del_z_neg: distance to extend grid in positive, negative x, y, and z directions from the center

		x, y, z_spacing: grid spacing in each dimension.


	3. Removal of ring piercings:
	
		num_rings: number of rings in the system.
		
		Ring_names_dict: specify the atom names for each ring manually
			e.g.: dict set Ring_names_dict 0 {C1 C2 C3 C4 C5 C6}

How to run:

	vmd 

	source run.tcl 
	
	Note: must first open vmd in gui mode and then source run.tcl

Output:

	dcds/dcdsstart_*.dcd: dcd containing each pose
