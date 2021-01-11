#!/bin/bash
for i in {1..8}
do 
~/Desktop/run-gpu/rungpu minimize_$i.namd 0 #### running each minimization namd script in parallel by submitting to multiple processors
done
