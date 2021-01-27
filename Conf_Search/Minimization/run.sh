#!/bin/bash
num_dcds=1 ## set number of dcds in system

for i in {1..$num_dcds}
do 
~/Desktop/run-gpu/rungpu minimize_$i.namd 0 #### running each minimization namd script in parallel by submitting to multiple processors
done
