******************************************************************************
Monte Carlo based Pathway Search (MCPS) for exploring molecular processes in high-dimensional conformational space
===============================================================================

Contributors: Archit Vasan and Nandan Haloi

Affiliation: University of Illinois at Urbana-Champaign

Emails: akvasan2@illinois.edu, nhaloi2@illinois.edu

About: 

Method to efficiently and systematically sample high-dimensional molecular processes while considering multiple slow degrees of freedom. In our implementation we focus on obtaining permeation pathways of antibiotics through outer membrane porins. 

Involves 3 steps:

1.  Exhaustive search for antibiotic poses within translational and rotational space.  After this search, a multidimensional energy landscape is created by evaluating antibiotic-protein interaction energy for each pose. (Found in Conf_Search directory) 

2.  Monte Carlo Based Pathway Search Algorith to walk through the energy landscape using Monte Carlo (MC) moves.  Since rotation and translation are slow degrees of freedom, limited changes in antibiotic orientation and position are allowed in each MC move. We advise that you run this multiple times to obtain multiple trajectories to sufficiently sample conformational space.  You can observe the convergence by plotting the trajectory density.  This code can be run on multiple  processors. (Found in MCPS directory)

3. Determination of most likely pathways sampled in our MCPS trajectories.  The trajectory data is used to construct a transition matrix which is inputted into Dijkstra's algorithm to obtain the most likely path. Can also be used to distinguish diverging paths. (Found in MostLikelyPathway directory)

The README files within each step of the method go into more detail about the necessary parameters.

Necessary softwares/programming environments:

	VMD
		Orient plugin  Instructions to install are at https://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/
	
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
	
