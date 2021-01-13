#!/bin/bash
for i in {1..8}
do 
namd2 +p3 pie_$i.namd > pie_$i.log &
done
