import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits import mplot3d
import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

##### Load in transition file ######
transition_file = '../MCPS/Output_Files/transition_search.dat' #### File containing transitions

t_file=open(transition_file)
transitions_tot=[]
with t_file as my_file:
    for line in my_file:
        myarray=np.fromstring(line, dtype=float, sep=' ')
        transitions_tot.append(myarray)
transitions=transitions_tot[0:2000]

data = np.loadtxt('../MCPS/Input_Files/input.dat') #### Input file from exhaustive searhc protocol

#### Grid files from grid construction step #####
z_grid = np.loadtxt('Grid_Data/grid_z.dat')
inc_grid = np.loadtxt('Grid_Data/grid_inc.dat')
az_grid = np.loadtxt('Grid_Data/grid_az.dat')

dens_grid = np.loadtxt('Grid_Data/grid_density.dat')
num_dens_grid = dens_grid

top_z = 4.5 #### Top and bottom of the grid
bott_z = -2.5

####### grids in z, inc, az spaces #########

min_z_grid = np.unique(z_grid)
min_inc_grid = np.unique(inc_grid)
min_az_grid = np.unique(az_grid)

###### Plot histogram of all densities for each grid #####
##### Helps in determining cutoff for density 

plt.hist(dens_grid)
plt.savefig('Images/dens_hist.png')
plt.close()
cmap="RdBu_r"

######## 3D plot of each density grid #######

fig=plt.figure()
ax = plt.axes(projection='3d')
dens = ax.scatter3D(z_grid,az_grid,inc_grid,c=num_dens_grid, cmap=cmap)

ax.set_xlim(bott_z,top_z)
ax.set_ylim(0,360)
ax.set_zlim(0,180)
plt.savefig('Images/gridding_density_all.png')
plt.close()

######### Plotting grids with density > cutoff ######

av_inc_z_dens = []
inc_grid_proj = []
z_grid_proj = []

dens_cut = 0.05 

z_grid_above_0p1 = z_grid[np.where(dens_grid>dens_cut)]
inc_grid_above_0p1 = inc_grid[np.where(dens_grid>dens_cut)]
az_grid_above_0p1 = az_grid[np.where(dens_grid>dens_cut)]


####### Clusters have azimuthal cutoff at 180 degrees

z_grid_c1 = z_grid_above_0p1[np.where(az_grid_above_0p1<180)]
inc_grid_c1 = inc_grid_above_0p1[np.where(az_grid_above_0p1<180)]
az_grid_c1 = az_grid_above_0p1[np.where(az_grid_above_0p1<180)]

z_grid_c2 = z_grid_above_0p1[np.where(az_grid_above_0p1>180)]
inc_grid_c2 = inc_grid_above_0p1[np.where(az_grid_above_0p1>180)]
az_grid_c2 = az_grid_above_0p1[np.where(az_grid_above_0p1>180)]

ax = plt.axes(projection='3d')
dens_c1 = ax.scatter3D(z_grid_c1,az_grid_c1,inc_grid_c1,c='black')#, cmap=cmap)
dens_c2 = ax.scatter3D(z_grid_c2,az_grid_c2,inc_grid_c2,c='forestgreen')#, cmap=cmap)

ax.set_xlim(-2,4.5)
ax.set_ylim(0,360)
ax.set_zlim(0,180)
plt.savefig('Images/gridding_density_clusters.png')
plt.close()

###### Plotting grids in z/inc spaces and z/az spaces ########

for z in range(len(min_z_grid)):
    for inc in range(len(min_inc_grid)):
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z])
        inc_grid_proj = np.append(inc_grid_proj,min_inc_grid[inc])
        bin_data = []
        for d in range(len(dens_grid)):
            if z_grid[d] == min_z_grid [z]:
                if inc_grid[d] == min_inc_grid[inc]:
                    bin_data = np.append(bin_data,dens_grid[d])
        av_inc_z_dens = np.append(av_inc_z_dens,np.mean(bin_data))

        z_grid_proj = np.append(z_grid_proj,min_z_grid[z]-0.25)     
        inc_grid_proj = np.append(inc_grid_proj,min_inc_grid[inc])
        av_inc_z_dens = np.append(av_inc_z_dens,np.mean(bin_data))
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z]-0.5)
        inc_grid_proj = np.append(inc_grid_proj,min_inc_grid[inc])
        av_inc_z_dens = np.append(av_inc_z_dens,np.mean(bin_data))
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z]-0.75)     
        inc_grid_proj = np.append(inc_grid_proj,min_inc_grid[inc])
        av_inc_z_dens = np.append(av_inc_z_dens,np.mean(bin_data))

plt.scatter(z_grid_proj,inc_grid_proj,c=av_inc_z_dens,cmap=cmap,s=500,marker='s')
plt.xlim(-2,4.5)
plt.ylim(0,180)
plt.savefig('Images/z_inc_proj.png')
plt.close()

av_az_z_dens = []
az_grid_proj = []
z_grid_proj = []

for z in range(len(min_z_grid)):
    for az in range(len(min_az_grid)):
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z])
        az_grid_proj = np.append(az_grid_proj,min_az_grid[az])
        bin_data = []
        for d in range(len(dens_grid)):
            if z_grid[d] == min_z_grid [z]:
                if az_grid[d] == min_az_grid[az]:
                    bin_data = np.append(bin_data,dens_grid[d])
        av_az_z_dens = np.append(av_az_z_dens,np.mean(bin_data))
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z]-0.25)     
        az_grid_proj = np.append(az_grid_proj,min_az_grid[az])
        av_az_z_dens = np.append(av_az_z_dens,np.mean(bin_data))
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z]-0.5)
        az_grid_proj = np.append(az_grid_proj,min_az_grid[az])
        av_az_z_dens = np.append(av_az_z_dens,np.mean(bin_data))
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z]-0.75)     
        az_grid_proj = np.append(az_grid_proj,min_az_grid[az])
        av_az_z_dens = np.append(av_az_z_dens,np.mean(bin_data))


plt.scatter(z_grid_proj,az_grid_proj,c=av_az_z_dens,cmap=cmap,s=500,marker='s')
plt.xlim(-2,4.5)
plt.ylim(0,360)
plt.savefig('Images/z_az_proj.png')
plt.close()


