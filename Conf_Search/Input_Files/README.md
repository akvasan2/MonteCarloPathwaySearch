**********************************************
Input Files 
===============================================================================
About:
	
	Input files used in Initial_Search, minimization, pair interaction energy, slow coordinate analysis
	make_cnst.tcl will make fix.pdb and PIE.pdb files which are important in minimization and pair interaction calculation 

Required files:
	
	Structure and coordinate (e.g. psf and pdb) files stored as system.psf and system.pdb
	
Important Parameters to specify in make_cnst.tcl: 
	
	sel1, sel2: 2 sets of atoms used when calculating pair interaction energy
	
	sel: atoms to be fixed during minimization

How to run:
	
	vmd make_cnst.tcl

Output:
	
	fix.pdb

	PIE.pdb
