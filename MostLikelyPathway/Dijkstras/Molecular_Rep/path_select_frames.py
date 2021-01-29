from __future__ import print_function, division, absolute_import
import csv
import numpy as np
import copy
import matplotlib.pyplot as plt
__all__ = ['paths', 'top_path']

import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

####### Load in data ###########3
###################################

t_file=open('../../../../MCPS/Output_Files/transition_search.dat')

transitions_tot=[]
with t_file as my_file:
    for line in my_file:
        myarray=np.fromstring(line, dtype=float, sep=' ')
        transitions_tot.append(myarray)
transitions=transitions_tot

data=np.loadtxt('../../../../MCPS/Input_Files/input.dat')

Most_likely_paths = np.loadtxt('../Cluster1/pathway.dat') ## pathway file obtained for Cluster 1 in the previous step.

top_z=np.max(data[:,0])
bott_z=np.min(data[:,0])

inc_step = 18
az_step = 18

z_step=1
z_num=int((top_z-bott_z)/float(1))

z_ind = 0
inc_ind = 1
az_ind = 2

output_file = "cluster1_path_frames.dat"

###### Create grids #############
################################

it=0

Grid_z=[]
Grid_ang=[]
Grid_az=[]

iteration=0

grid_z=[]

for i in range(z_num):
    z_i=top_z-z_step*i
    grid_z=np.append(grid_z,z_i)

grid_inc_ang=np.arange(0,180,inc_step)
grid_az_ang=np.arange(0,360,az_step)

for z in range(len(grid_z)):
    for ang in range(len(grid_inc_ang)):
        for az in range(len(grid_az_ang)): 
            Grid_z=np.append(Grid_z,grid_z[z])
            Grid_ang=np.append(Grid_ang,grid_inc_ang[ang])
            Grid_az=np.append(Grid_az,grid_az_ang[az])
            iteration+=1

colors=['black','gray','orange','red','green']

################ Assign Grids to most likely pathway ############
#################################################################

mlp_z=[]
mlp_ang=[]
mlp_az=[]

end_path = Most_likely_paths[-1]

for j in range(len(Most_likely_paths)):
    if Most_likely_paths[j]!=end_path and Most_likely_paths[j]!=0:
        iteration_number=Most_likely_paths[j] - 2
        mlp_z=np.append(mlp_z,Grid_z[int(iteration_number)])
        mlp_ang=np.append(mlp_ang,Grid_ang[int(iteration_number)])
        mlp_az=np.append(mlp_az,Grid_az[int(iteration_number)])
MLP_z=mlp_z
MLP_ang=mlp_ang
MLP_az=mlp_az
Most_likely_path_frames_comp={}
MLP_low_en_frame=[]

############### ID frames in each bin of path ##############
#################################################################3

f=open(output_file,'w')
f.write("")
f.close()

for i in range(0,len(MLP_z)):
    Most_likely_path_frames=[]
    Most_likely_path_frames_comp[i]=[]
    mlp_z=MLP_z[i]
    mlp_ang=MLP_ang[i]
    mlp_az=MLP_az[i]
    for p in range(len(transitions)):
        for f in range(len(transitions[p])):
            ind=int(transitions[p][f])
            if data[ind][z_ind]<=mlp_z and data[ind][z_ind]>=mlp_z-1 and data[ind][inc_ind]>=mlp_ang and data[ind][inc_ind]<=mlp_ang+inc_step and data[ind][az_ind]>=mlp_az and data[ind][az_ind]<=mlp_az+az_step:
                Most_likely_path_frames_comp[i]=np.append(Most_likely_path_frames_comp[i],ind)
                Most_likely_path_frames=np.append(Most_likely_path_frames,ind)
    f=open(output_file,'a+',newline='')
    writer=csv.writer(f,delimiter=' ',quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(Most_likely_path_frames_comp[i])
    f.close()

