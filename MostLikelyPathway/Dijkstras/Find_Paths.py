import numpy as np
import math as m
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import BuildTransitionMatrix as BTM
import dijkstra
######################### Opening transition file with filtered cluster1 paths##################################

t_file=open('../Traj_Group_Data/cluster1_paths.dat')
transitions_tot=[]
with t_file as my_file:
    for line in my_file:
        myarray=np.fromstring(line, dtype=float, sep=' ')
        transitions_tot.append(myarray)
transitions_c1=transitions_tot

######################### Opening transition file with filtered cluster2 paths##################################

t_file=open('../Traj_Group_Data/cluster2_paths.dat')
transitions_tot=[]
with t_file as my_file:
    for line in my_file:
        myarray=np.fromstring(line, dtype=float, sep=' ')
        transitions_tot.append(myarray)

transitions_c2=transitions_tot


######################## Input data from docking #######################3
data=np.loadtxt('../../Input_Files/input.dat')

###################Defining parameters #########################
####################D###################D##################D

z_ind = 0
inc_ind = 1
az_ind = 2
en_ind = 3 

beta=np.std(data[:,en_ind])

##### Define z bins

top_z=np.max(data[:,0])
print(top_z)
bott_z=-11.5
print(bott_z)

##### Step in each z, inc, az bin

z_step=1 
inc_step = 18
az_step = 18

z_num=int((top_z-bott_z)/float(z_step))
bins_z=[]

for i in range(z_num):
    z_i=top_z-z_step*i
    bins_z=np.append(bins_z,z_i)

##### Define angular bins

bins_inc_ang=np.arange(0,180,inc_step)
bins_az_ang=np.arange(0,360,az_step)

Met_mat_file_c1 = 'Cluster1/Data_Files/met_matrix_grid.dat' #### where to store cluster1 metropolis matrix

Neg_log_file_c1 = 'Cluster1/Data_Files/neg_log_matrix_grid.dat' ##### where to store cluster 1 neg log probability matrix

Met_mat_file_c2 = 'Cluster2/Data_Files/met_matrix_grid.dat' #### where to store cluster1 metropolis matrix

Neg_log_file_c2 = 'Cluster2/Data_Files/neg_log_matrix_grid.dat' ##### where to store cluster 1 neg log probability matrix

######################### Building transition matrix ########################

#### cluster 1
BTM.BuildTransitionMatrix(transitions_c1, data, z_ind, inc_ind, az_ind, en_ind, beta, top_z, bott_z, z_step, inc_step, az_step, z_num, bins_z, bins_inc_ang, bins_az_ang, Met_mat_file_c1, Neg_log_file_c1)

neg_log_matrix_c1=np.loadtxt('Cluster1/Data_Files/neg_log_matrix_grid.dat')

graph_c1={index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(neg_log_matrix_c1)}

start=0
end=len(neg_log_matrix_c1)-1

short_dist, track_path = dijkstra.dijkstra(graph_c1, start,end)

np.savetxt('Cluster1/pathway.dat',track_path)

#### cluster 2

BTM.BuildTransitionMatrix(transitions_c2, data, z_ind, inc_ind, az_ind, en_ind, beta, top_z, bott_z, z_step, inc_step, az_step, z_num, bins_z, bins_inc_ang, bins_az_ang, Met_mat_file_c2, Neg_log_file_c2)

neg_log_matrix_c2=np.loadtxt('Cluster2/Data_Files/neg_log_matrix_grid.dat')

graph_c2={index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(neg_log_matrix_c2)}

start=0
end=len(neg_log_matrix_c2)-1

short_dist, track_path = dijkstra.dijkstra(graph_c2, start,end)

np.savetxt('Cluster2/pathway.dat',track_path)
