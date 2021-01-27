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

########## Load in transition data ############

t_file=open('../../../Output_Files/transition_search.dat')
transitions_tot=[]
with t_file as my_file:
    for line in my_file:
        myarray=np.fromstring(line, dtype=float, sep=' ')
        transitions_tot.append(myarray)
transitions=transitions_tot

z_step = 1
inc_step = 18
az_step = 18

data=np.loadtxt('../../../Input_Files/input.dat')

###### Load in most likely path from Dijkstra's algorithm ######

Most_likely_paths=[int(p) for p in np.loadtxt('pathway.dat')]

###### Create grids in z, incliantion and azimuthal spaces.  These grids are used to plot the eventual path ######

it=0

Grid_z=[]
Grid_ang=[]
Grid_az=[]
iteration=0

top_z=np.max(data[:,0])
bott_z=np.min(data[:,0])
z_step=1
z_num=int((top_z-bott_z)/float(z_step))
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

######## For each point along the path, determine the associated grid value ############

mlp_z=[]
mlp_ang=[]
mlp_az=[]

for j in range(len(Most_likely_paths)):
    if Most_likely_paths[j]!=7402 and Most_likely_paths[j]!=0:
        iteration_number=Most_likely_paths[j]
        mlp_z=np.append(mlp_z,Grid_z[int(iteration_number)])
        mlp_ang=np.append(mlp_ang,Grid_ang[int(iteration_number)]+9)
        mlp_az=np.append(mlp_az,Grid_az[int(iteration_number)]+9)

MLP_z=mlp_z
MLP_ang=mlp_ang
MLP_az=mlp_az

####### Plot the paths !#########

plt.plot(MLP_z,MLP_ang,color='black')
plt.scatter(MLP_z,MLP_ang,color='black')
plt.ylim(0,180)
plt.savefig(f'path_scatter_inc_z_full.png',bbox_inches='tight',transparent=True)
plt.close()

plt.plot(MLP_z,MLP_az,color='black')
plt.scatter(MLP_z,MLP_az,color='black')
plt.ylim(0,360)
plt.savefig(f'path_scatter_az_z_full.png',bbox_inches='tight',transparent=True)
plt.close()

