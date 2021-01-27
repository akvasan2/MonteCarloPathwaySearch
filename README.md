******************************************************************************
Monte Carlo based Pathway Search (MCPS) for exploring molecular processes in high-dimensional conformational space
===============================================================================

Contributors: Archit Vasan and Nandan Haloi

Affiliation: University of Illinois at Urbana-Champaign

Emails: akvasan2@illinois.edu, nhaloi2@illinois.edu

About: 

Method to efficiently and systematically sample high-dimensional molecular processes while considering multiple slow degrees of freedom. In our implementation, we focus on obtaining permeation pathways of antibiotics through outer membrane porins. 

Involves 3 steps:

1.Exhaustive search for antibiotic poses within translational and rotational space.  After this search, a multidimensional energy landscape is created by evaluating antibiotic-protein interaction energy for each pose. (Found in Conf_Search directory) 

   Divided into 4 substeps:

A.Exhaustive initial search of all poses for the antibiotic. (Found in Conf_Search/Initial_Search)
Involves:
	
	Generating multiple drug orientations by aligning a vector of the drug to fibonacci sphere and doing self rotations of the drug
	
	Translating the drug to all possible positions within the protein at points along a grid.

	Removing any clashes or ring pierces between drug and protein.
		

Required files:
        	
	Structure and coordinate files (i.e. psf, pdb) files stored as ../Input_Files/system.psf and ../Input_Files/system.pdb

How to run:
		
	vmd

	source run.tcl

	Note: must first open vmd in gui mode and then source run.tcl

Output:
	dcds/dcdsstart_*.dcd: dcd containing each pose

B.Minimization of each pose obtained from initial search.  (Found in Conf_Search/Minimization) 
	
Required files:

	Structure and coordinate files (i.e. psf, pdb) files stored as ../Input_Files/system.psf and ../Input_Files/system.pdb

	Parameter files: stored in ../../Parameters
        
        fixedAtomsFile specifying atoms to be fixed: stored as ../../Input_Files/fix.pdb. The B column of each atom of the two molecules should equal 1 or 2 
        
How to run:

	bash run.sh

Output:

	minimize_output/dcdsstart_*.dcd


C. Evaluation of pair interaction energy of the drug poses. (Found in Analysis/PIE)

D. Evaluation of slow coordinates of the drug poses.  (Found in Analysis/Slow_Coordinates) 


2. Monte Carlo Based Pathway Search (MCPS) Algorithm to walk through the energy landscape using MC moves. Since rotation and translation are slow degrees of freedom, limited changes in antibiotic orientation and position are allowed in each MC move. Need to run this multiple times to obtain multiple trajectories such that interested conformational space is sufficiently sampled. You can determine the convergence by plotting the trajectory density, projected onto the individual conformation spaces. This code can be run on multiple processors. (Found in MCPS directory) 

3. Determination of most likely pathways sampled in our MCPS trajectories. The trajectory data is used to construct a transition matrix which is inputted into Dijkstra's algorithm to obtain the most likely path. Can also be used to distinguish diverging paths. (Found in MostLikelyPathway directory)

The README files within each step of the method go into more details about the necessary parameters.

Necessary softwares/programming environments:

	VMD
	Additional plugins necessary: 
		Orient plugin Instructions to install are at https://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/
	
	Python 3
	Modules necessary:
		numpy
		math
		random
		multiprocessing
		joblib
		csv
		matplotlib
	NAMD2
