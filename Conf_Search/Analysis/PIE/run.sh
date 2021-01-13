#!/bin/bash
num_dcds=1 ## set number of dcds in system

#### create config files for each dcd
for i in $(seq 2 1 $num_dcds)
do
       
        cp pie_1.namd pie_${i}.namd
        sed -i "s/start_1/start_${i}/g" pie_${i}.namd

done

#### run calculation for each dcd
for i in $(seq 1 1 $num_dcds)
do 
	#### store each log file in Log_Files/Raw
	namd2 +p3 pie_$i.namd > Log_Files/Raw/pie_$i.log &
done

#### clean Log_Files directory
rm -f Log_Files/*.log

#### extract important energy information from each log file
for i in $(seq 1 1 $num_dcds)

do
	ipython extract.py Log_Files/Raw/pie_${i}.log Log_Files/pie_${i}.log
done

### merge each log file into PIE.dat
rm -f PIE.dat

cat Log_Files/*.log >> Output/PIE.dat
