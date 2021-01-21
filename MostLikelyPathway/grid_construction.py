import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)
cmap="RdBu_r"

########### Define functions ########

def grid_construction (data, transitions, z_index, inc_index, az_index, en_index, en_mean, en_std, top_z, bott_z, z_step, inc_step, az_step) :

	 ###### Store transition z, inclination, azimuthal, boltzmann weights if z for data is within grid values (within bott_z and top_z).

	transitions_angles=[]
	transitions_z=[]
	transitions_beta=[]
	transitions_weights=[]
	for i in range(len(transitions)):
		for j in range(len(transitions[i])-1):
			f=int(transitions[i][j])
			#### store data from trajectory into a large array if it lies within bott_z and top_z
			if data[f][z_index]>bott_z and data[f][z_index]<top_z:
				transitions_z=np.append(transitions_z,data[f][z_index])
				transitions_angles=np.append(transitions_angles,data[f] [inc_index])
				transitions_beta=np.append(transitions_beta,data[f][az_index])
				transitions_weights=np.append(transitions_weights,np.exp(-(data[f][en_index]-en_mean)/float(en_std)))
	
	####### Create a list of all z, incliantion, azimuthal bins

	bins_z=[]

	z_num=int((top_z-bott_z)/float(z_step))

	for i in range(z_num):
	    z_i=top_z-z_step*i
	    bins_z=np.append(bins_z,z_i)
	bins_inc_ang=np.arange(0,180,inc_step)#### Bins for inclination, azimuthal such that number specifies lower limit of grid cell.
	bins_az_ang=np.arange(0,360,az_step)

	###### This next step iterates over bins to create a 3d grid

	##### These grid arrays essentially are a 3D grid split up into the corresponding columns
	Grid_z = []
	Grid_inc = []
	Grid_az = []
	Len_Grid = (len(bins_z))*len(bins_inc_ang)*len(bins_az_ang) #### Length of the full 3D grid
		
	Grid_density = np.zeros(Len_Grid) #### Density of each point in grid.  Initialized to zero at each point right now
	it = 0

	for z in range(len(bins_z)):
	    for inc in range(len(bins_inc_ang)):
	        for az in range(len(bins_az_ang)):
	            it+=1
	            Grid_z=np.append(Grid_z,bins_z[z])
	            Grid_inc=np.append(Grid_inc,bins_inc_ang[inc])
	            Grid_az=np.append(Grid_az,bins_az_ang[az])
	            for n in range(len(transitions_z)): #### all stored data from trajectories.
			#### check if stored transitions_z, transitions_inc, transitions_az at n are within the correct bin
	                if transitions_z[n]<bins_z[z] and transitions_z[n]>bins_z[z]-z_step:
	                    if transitions_angles[n]>bins_inc_ang[inc] and transitions_angles[n]<bins_inc_ang[inc]+inc_step:
	                        if transitions_beta[n]>bins_az_ang[az] and transitions_beta[n]<bins_az_ang[az]+az_step:
	                            Grid_density[it-1]+=transitions_weights[n] # based on stored value in transitions_weights array, increase density at this grid.
	
	return (Grid_z,Grid_inc,Grid_az,[i/float(z_step*inc_step*az_step) for i in Grid_density]) ### Return Grid data for each dimension.  normalize density by the size of each dimension 


