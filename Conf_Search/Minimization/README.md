===============================================================================
Minimization of each pose from exhaustive sampling procedure
===============================================================================

About:
	Scripts used to minimize all poses created from exhaustive search in ../Initial_Search 
Required input files:
	structure and coordinate (i.e. psf and pdb) files: stored as ../../Input_Files/system.psf and ../../Input_Files/system.pdb
	parameter files: stored in ../../Parameters
	fixedAtomsFile specifying atoms to be fixed: stored as ../../Input_Files/fix.pdb
		B column of each atom of the two molecules should equal 1 or 2 
	initially placed pose dcd files: stored as ../Initial_Search/dcds/dcdsstart_*.dcd

Important Parameters to specify in minimize_1.namd: 
	structure and coordinate files names 
		(e.g. ../../Input_Files/system.psf and ../../Input_Files/system.pdb)
	whether to use GBIS minimization or not
	minimization time

Important Parameters to specify in run.sh: 
	num_dcds: number of dcd files

How to run:
	bash run.sh 

Output:
	minimize_output/dcdsstart_*.dcd
