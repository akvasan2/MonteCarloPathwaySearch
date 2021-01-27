import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv

import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

cmap="RdBu_r"

######### Functions Defined Below #########
########################################################################

########################################################################

######### For each trajectory, assign each point to a particular grid ############

def assign_paths(transitions, data, z_ind, inc_ind, az_ind, top_z,bott_z,bins_z,bins_inc_ang,bins_az_ang,paths_file,paths_density_file,z_step,inc_step, az_step,Grid_z,Grid_density):
        
    pathway_grids = {}
    
    pathway_density = {}
    
    pathway_z = {}
        
    for i in range(len(transitions)):
            pathway_grids[i] = []
            pathway_density[i] = []
            for j in range(len(transitions[i])):
                f = int(transitions[i][j])
                it = 0
                #### Assign point along pathway to a particular grid 
                for z in range(len(bins_z)):
                    if data[f][z_ind]<=bins_z[z] and data[f][z_ind]>=bins_z[z]-z_step:
                        for inc in range(len(bins_inc_ang)):
                            if data[f][inc_ind]>=bins_inc_ang[inc] and data[f][inc_ind]<=bins_inc_ang[inc]+inc_step:
                                for az in range(len(bins_az_ang)):
                                    it = (z*len(bins_inc_ang)*len(bins_az_ang)) + inc*len(bins_az_ang) + az+1 
                                    if data[f][az_ind]>=bins_az_ang[az] and data[f][az_ind]<=bins_az_ang[az]+az_step:
                                        pathway_grids[i] = np.append( pathway_grids[i], it)
                                        print(it)
                                        print(Grid_density[it-1])
                                        pathway_density[i] = np.append( pathway_density[i], Grid_density[it-1] )
                                        break

    for i in range(len(pathway_grids)):
        
        f=open(paths_file,'a+',newline='')
        writer=csv.writer(f,delimiter=' ',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(pathway_grids[i])
        f.close()                                                                          
        
        f=open(paths_density_file,'a+',newline='')
        writer=csv.writer(f,delimiter=' ',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(pathway_density[i])
        f.close()                                                                          
    return (pathway_grids,pathway_density)    
    
####### Use if you have already assigned paths to grids

def open_paths (paths_file,paths_density_file):
    p_file=open(paths_file)
    pathway_grids=[]
    with p_file as my_file:
        for line in my_file:
            myarray=np.fromstring(line, dtype=float, sep=' ')
            pathway_grids.append(myarray)
    
    p_file=open(paths_density_file)
    pathway_density=[]
    with p_file as my_file:
        for line in my_file:
            myarray=np.fromstring(line, dtype=float, sep=' ')
            pathway_density.append(myarray)

    return (pathway_grids,pathway_density)

#################################################################################################################

def cluster_trajectories(transitions, bins_z, pathway_grids, pathway_density, count_cutoff, density_cutoff, cluster_cutoff,paths_file,paths_file_c1,paths_file_c2, path_gridfle_c1, path_gridfle_c2,Grid_z,Grid_az,Grid_density):

    allowed_trajectories = {}
    allowed_path_indices = []
    it = -1
    it_c1 = -1
    it_c2 = -1
    clust1_traj = {}
    clust2_traj = {}
    clust1_grid_traj = {}
    clust2_grid_traj = {}

    for i in range(len(transitions)):
        count = 0  
        clust1_cnt = 0  #### number of points along the trajectory that belong to cluster 1.  set to 0 initially
        clust2_cnt = 0 #### number of points along the trajectory that belong to cluster 2.  set to 0 initially

        #### In order to group trajectories, check that at every z, at least one point of the pathway lies within that group and not in any other defined groups

        for z in bins_z:
            countz = 0
            clust1_cnt_z = 0 #### Initialize this to 0 but if a point along the trajectory lies within clust1, then increase this value 
            clust2_cnt_z = 0 #### Do the same as the above comment but for clust2 

            ##### Run through the points of each trajectory (which have been assigned to grids already) and determine whether that point belongs to clust1 or clust2

            for j in range(len(pathway_grids[i])):
                f = int(pathway_grids[i][j])
                if Grid_z[f-1]==z and pathway_density[i][j]>=density_cutoff: #### Only use this point if the density is > density_cutoff and the z lies within the correct location 
                    countz += 1
                    ###### Depending on if the point belongs to a particular cutoff, increase the counter 
                    if Grid_az[f]<cluster_cutoff:
                        clust1_cnt_z+=1
                    else:
                        clust2_cnt_z+=1
            if countz !=0:
                count+=1
                if clust1_cnt_z !=0 and clust2_cnt_z == 0: # If at this z, all classified points are cluster1 and not clust2, increase this count
                    clust1_cnt+=1
                elif clust2_cnt_z !=0 and clust1_cnt_z == 0:
                    clust2_cnt+=1
        ##### If at every z, at least one trajectory point is within a particular cluster then cluster the trajectory  
        if count>=count_cutoff:
            it+=1
            allowed_trajectories[it] = transitions[i]
            if clust1_cnt>=count_cutoff:
                it_c1 +=1
                clust1_traj[it_c1] = transitions[i]
                clust1_grid_traj[it_c1] = pathway_grids[i] 

            elif clust2_cnt>=count_cutoff:
                it_c2 +=1
                clust2_traj[it_c2] = transitions[i]
                clust2_grid_traj[it_c2] = pathway_grids[i] 
###############################################################

##### Saving allowed trajectories

##############################################
        
    f = open(paths_file, 'w')
    f.write('')
    f.close()
    
    for i in range(len(allowed_trajectories)):
        
        f=open(paths_file,'a+',newline='')
        writer=csv.writer(f,delimiter=' ',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(allowed_trajectories[i])
        f.close()                                                                          
    
    f = open(paths_file_c1,'w')
    f.write('')
    f.close()
    
    for i in range(len(clust1_traj)):
        
        f=open(paths_file_c1,'a+',newline='')
        writer=csv.writer(f,delimiter=' ',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(clust1_traj[i])
        f.close()                                                                          
    
    f = open(paths_file_c2,'w')
    f.write('')
    f.close()
    
    for i in range(len(clust2_traj)):
        
        f=open(paths_file_c2,'a+',newline='')
        writer=csv.writer(f,delimiter=' ',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(clust2_traj[i])
        f.close()                                                                          

    f = open(path_gridfle_c1,'w')
    f.write('')
    f.close()
    
    for i in range(len(clust1_grid_traj)):
        
        f=open(path_gridfle_c1,'a+',newline='')
        writer=csv.writer(f,delimiter=' ',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(clust1_grid_traj[i])
        f.close()                                                                          

    f = open(path_gridfle_c2,'w')
    f.write('')
    f.close()
    
    for i in range(len(clust2_grid_traj)):
        
        f=open(path_gridfle_c2,'a+',newline='')
        writer=csv.writer(f,delimiter=' ',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(clust2_grid_traj[i])
        f.close()                                                                          



