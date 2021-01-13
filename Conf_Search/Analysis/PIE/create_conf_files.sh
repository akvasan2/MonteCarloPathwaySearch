######## num_dcds:how many dcd files stored in Initial_Search 
num_dcds=8
for i in $(seq 2 1 $num_dcds)
do
       
        cp pie_1.namd pie_${i}.namd
        sed -i "s/start_1/start_${i}/g" pie_${i}.namd

done

