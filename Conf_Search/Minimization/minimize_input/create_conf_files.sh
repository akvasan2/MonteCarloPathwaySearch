######## num_dcds:how many dcd files stored in Initial_Search 
num_dcds=8
for i in $(seq 2 1 $num_dcds)
do
       
        cp minimize_1.namd minimize_${i}.namd
        sed -i "s/start_1/start_${i}/g" minimize_${i}.namd

done

