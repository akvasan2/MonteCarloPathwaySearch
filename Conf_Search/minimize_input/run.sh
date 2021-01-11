#!/bin/bash
for i in {1..8}
do 
~/Desktop/run-gpu/rungpu minimize_$i.namd 0
done
