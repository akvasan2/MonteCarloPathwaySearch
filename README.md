******************************************************************************
# Monte Carlo based Pathway Search (MCPS) for exploring molecular processes in high-dimensional conformational space

Contributors: Archit Vasan and Nandan Haloi

Affiliation: University of Illinois at Urbana-Champaign

Emails: akvasan2@illinois.edu, nhaloi2@illinois.edu

If you want to use MCPS please cite: Chem. Sci., 2021,12, 15028-15044
## General overview: 

Method to efficiently and systematically sample high-dimensional molecular processes while considering multiple slow degrees of freedom. In our implementation, we focus on obtaining permeation pathways of antibiotics through outer membrane porins. 

## Protocol involves 3 steps:

### 1. Exhaustive search for antibiotic poses (Found in Conf_Search directory) 


#### Code description:

Exhaustively searches for antibiotic poses within translational and rotational space.  After this search, a multidimensional energy landscape is created by evaluating antibiotic-protein interaction energy for each pose. 

#### Divided into 4 substeps:

#### A. Exhaustive initial search of all poses for the antibiotic. 

Involves: generating multiple drug orientations, translating the drug to all possible positions within the protein at points along a grid, and removing clashes/ring pierces between drug and protein. (Found in Conf_Search/Initial_Search) 

##### Required files:
        	
	Structure and coordinate files (i.e. psf, pdb) files stored as Conf_Search/Input_Files/system.psf and Conf_Search/Input_Files/system.pdb

##### How to run:
	
	Set parameters in run.tcl 
	
	vmd

	source run.tcl

	Note: must first open vmd in gui mode and then source run.tcl
	
##### Output:

	dcds/dcdsstart_*.dcd: dcd containing each pose

#### B. Minimization of each pose obtained from initial search.  (Found in Conf_Search/Minimization) 
	
##### Required files:

	Structure and coordinate files (i.e. psf, pdb) files stored as Conf_Search/Input_Files/system.psf and Conf_Search/Input_Files/system.pdb

	Parameter files: stored in Conf_Search/Parameters
        
	fixedAtomsFile specifying atoms to be fixed: stored as ../../Input_Files/fix.pdb. The B column of each atom of the two molecules should equal 1 or 2 
        
	fix.pdb can be created by running: vmd make_cnst.tcl in the Conf_Search/Input directory.

##### How to run:

	Set parameters in minimize_1.namd and run.sh
	
	bash run.sh
	
##### Output:

	minimize_output/dcdsstart_*.dcd

#### C. Evaluation of pair interaction energy of the minimized drug poses. (Found in Conf_Search/Analysis/PIE)

##### Required files:
	
	Structure and coordinate (e.g. psf and pdb) files stored in Conf_Search/Input_Files
	
	Parameter files: stored in Conf_Search/Parameters
	
	pairInteractionFile (PIE.pdb):specifying atoms of the pair of molecules for calculating PIE.  B column of each atom of the two molecules should equal 1 or 2 
	PIE.pdb can be created by running: vmd make_cnst.tcl in the Conf_Search/Input directory.

##### How to run:
	
	Set parameters in pie_1.namd and run.sh
	
	bash run.sh 
		
##### Output:
	
	Output/PIE.dat

#### D. Evaluation of values for the slow coordinates (z-coordinate, inclination, azimuthal) for the drug poses. (Found in Conf_Search/Analysis/Slow_Coordinates) 

##### Required files:
	
	Structure (i.e. psf) file stored as Conf_Search/Input_Files/system.psf

##### How to run:
	
	Set parameters in slow_coor.tcl
	
	vmd slow_coor.tcl 

##### Output:

	slow_coor.dat: columns are z-coor, inclination, azimuthal data for each pose 

### 2. Monte Carlo Based Pathway Search (MCPS). (Found in MCPS directory) 

#### Code description:

Algorithm to walk through the derived multi-dimensional energy landscape using MC moves. Since rotation and translation are slow degrees of freedom, limited changes in antibiotic orientation and position are allowed in each MC move. Need to run this multiple times to obtain multiple trajectories such that conformational space is sufficiently sampled. You can determine the convergence by plotting the trajectory density, projected onto the individual conformation spaces. This code can be run on multiple processors. 

##### How to run:

	Set parameters in run.py
	
	python run.py
	
##### Output:

	Output_Files/transition_search.dat


### 3. Most Likely Pathway Identification. (Found in MostLikelyPathway directory)


#### Code description:

Determination of most likely pathways sampled in our MCPS trajectories. The trajectory data is used to construct a transition matrix which is inputted into Dijkstra's algorithm to obtain the most likely path. Can also be used to distinguish diverging paths. In our implementation, we use this method to obtain two diverging pathways; however, it can also be used to obtain either one or greater than 2 paths, depending on specific cases.

#### Divided into 3 substeps:

#### A. Filtering trajectories (Found in MostLikelyPathway)	

##### How to run:

	Set parameters in Filter_Trajectories.py
	
	python Filter_Trajectories.py
	
##### Output:

	Traj_Group_Data/cluster1_paths.dat
	   
	Traj_Group_Data/cluster2_paths.dat

#### B. Identification of most likely pathways (Found in MostLikelyPathway/Dijkstras) 

##### How to run:

	Set parameters in Dijkstras/Find_Paths.py
	
	python Dijkstras/Find_Paths.py
	
##### Output:

	Cluster1/pathway.dat: Most Likely Path for each group of trajectories, specifying the grid at each point along the path  

	Cluster2/pathway.dat

##### Plot the paths (Cluster*/pathway.dat) using plotting_path.py

#### C. Determination of representative snapshots for each path

##### How to run:

	Set parameters in path_select_frames.py and make_trajectory.tcl
	
	python path_select_frames.py
	
	vmd make_trajectory.tcl

##### Output:
	path_frames.dat:  all possible frames at each point along the path
	accepted_pathway.dat : frames for points along the final pathway obtained in algorithm.
	PDBs/window_*.pdb: snapshot at each point along path.

## Necessary softwares/programming environments:

### VMD

#### Additional VMD plugins necessary: 

Orient: Instructions to install are at https://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/

### Python 3

#### Modules necessary:

numpy

math

random

multiprocessing

joblib

csv

matplotlib

### NAMD2
