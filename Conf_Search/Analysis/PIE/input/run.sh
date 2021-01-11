#!/bin/bash
for i in {1..8}
do 
namd2 +p3 pie_diel_$i.namd > pie_diel_$i.log &
done
