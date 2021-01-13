**********************************************
Exhaustive search for antibiotic poses in conformational space 
===============================================================================


Scripts used to exhaustively search over the translational and rotational space for an antibiotic.  Poses determined are minimized and their interaction energy and slow degrees of freedom are evaluated.  

This protocol generates data that is crucial to running the MCPS algorithm to determine most likely pathways.

The order of steps to use this protocol are:
	
	1. An exhaustive initial search of all poses for the antibiotic.  Scripts found in Initial_Search/
	2. Minimization of each pose.  Scripts found in Minimization/ 
	3. Evaluation of pair interaction energy. Scripts found in Analysis/PIE/
	4. Evaluation of slow coordinates of the drug.  Scripts found in Analysis/Slow_Coordinates/ 

Input_Files and Parameters directories include structures and parameters needed to run these scripts, respectively.

Follow instructions detailed in the README.MD file within the specified directory for each step . 
