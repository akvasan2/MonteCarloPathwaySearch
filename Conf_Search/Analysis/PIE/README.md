
Evaluate Pair Interaction Energy (PIE) for each minimized pose

About:
	
	Scripts used to analyze the pair interaction energy for each pose which was
	previously minimized

Required files:
	
	structure and coordinate (e.g. psf and pdb) files stored in ../../Input_Files
	
	parameter files: stored in ../../Parameters
	
	pairInteractionFile (PIE.pdb):specifying atoms of the pair of molecules calculating PIE for
	
		B column of each atom of the two molecules should equal 1 or 2 
	
	all minimized pose dcd files stored in ../../Minimization/minimize_output 

Important Parameters to specify in pie_1.namd: 
	
	structure and coordinate files 
		
		(e.g. ../../Input_Files/system.psf and ../../Input_Files/system.pdb)
	
	dielectric: dielectric constant to run simulations default: 78.5	

Important Parameters to specify in pie_1.namd: 
	
	num_dcds: number of dcd files

How to run:
	
	bash run.sh 

Output:
	
	Output/PIE.dat
