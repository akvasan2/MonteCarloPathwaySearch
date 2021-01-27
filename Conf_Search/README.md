**********************************************
Exhaustive search for antibiotic poses in conformational space 
===============================================================================

Scripts used to exhaustively search over the translational and rotational space for an antibiotic to generate multiple poses.  Identified poses are minimized and their interaction energy and slow degrees of freedom are evaluated.  

This protocol generates a conformational landscape which MCPS walks through using Monte Carlo moves.

The order of steps to use this protocol are:
	
1. An exhaustive initial search of all poses for the antibiotic.  (Found in Initial_Search)

2. Minimization of each pose.  (Found in Minimization) 

3. Evaluation of pair interaction energy of the drug poses. (Found in Analysis/PIE)

4. Evaluation of slow coordinates of the drug poses.  (Found in Analysis/Slow_Coordinates) 

Follow instructions detailed in the README.MD file within each step's specified directory. 

Input_Files and Parameters directories include structures and parameters needed to run scripts, respectively.
