import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib
import kde_2d_weighted as KDE_2D

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)
cmap="RdBu_r"

##############################################################################################################################
############################################loading in data#########################################################

t_file=open('Output_Files/transition_search.dat')

transitions_tot=[]
with t_file as my_file:
    for line in my_file:
        myarray=np.fromstring(line, dtype=float, sep=' ')
        transitions_tot.append(myarray)
transitions=transitions_tot[0:2000]

data=np.loadtxt('Input_Files/input.dat')
beta=1/float(np.std(data[:,3]))
en_mn=np.mean(data[:,3])

transitions_z={}
transitions_beta={}
transitions_angles={}

#############################################################################################################################
####################################################creating arrays#########################################################


z_ecutoff_0=[]
a1_ecutoff_0=[]
a2_ecutoff_0=[]
ecutoff_0=[]

accepted=0
transitions_angles=[]
transitions_z=[]
transitions_beta=[]
transitions_energy=[]

for i in range(len(transitions)):
    for j in range(len(transitions[i])-1):
        f=int(transitions[i][j])
        if data[f][0]>-7.5 and data[f][0]<16.5:
            transitions_angles=np.append(transitions_angles,data[f] [1])
            transitions_z=np.append(transitions_z,data[f][0])
            transitions_beta=np.append(transitions_beta,data[f][2])
            transitions_energy=np.append(transitions_energy,data[f][3])

cmap="RdBu_r"


##############################################################################################################################

################################  Creating grids #################################################3

grid_z=[]
top_z=16.5
bott_z=-7.5
z_step=1
z_num=int((top_z-bott_z)/float(z_step))

for i in range(z_num):
    z_i=top_z-z_step*i
    grid_z=np.append(grid_z,z_i)

inc_step = 18
az_step = 18

grid_inc_ang=np.arange(0,180,inc_step)
grid_az_ang=np.arange(0,360,az_step)

Grid_z = []
Grid_inc = []
Grid_az = []
iteration = 0

for z in range(len(grid_z)):
    for inc in range(len(grid_inc_ang)):
        for az in range(len(grid_az_ang)):
            iteration+=1

Grid_density = np.zeros(iteration)

it = 0
for z in range(len(grid_z)):
    for inc in range(len(grid_inc_ang)):
        for az in range(len(grid_az_ang)):
            it+=1
            Grid_z=np.append(Grid_z,grid_z[z])
            Grid_inc=np.append(Grid_inc,grid_inc_ang[inc])
            Grid_az=np.append(Grid_az,grid_az_ang[az])
            frames = np.where((transitions_z<grid_z[z]) & (transitions_z>grid_z[z]-z_step) & (transitions_angles>grid_inc_ang[inc]) & \
                    (transitions_angles<grid_inc_ang[inc]+inc_step) & (transitions_beta>grid_az_ang[az]) & \
                    (transitions_beta<grid_az_ang[az]+az_step))
            Grid_density[it-1]=np.sum([np.exp(-beta*(transitions_energy[n]-en_mn)) for n in frames])

np.savetxt('grid_density.dat',[i/float(inc_step*az_step) for i in Grid_density])
np.savetxt('grid_z.dat',Grid_z)
np.savetxt('grid_inc.dat',Grid_inc)
np.savetxt('grid_az.dat',Grid_az)

##############################################################################################################################

################################  Creating projected density #################################################3
min_z_grid = grid_z
min_inc_grid = grid_inc_ang
min_az_grid = grid_az_ang
dens_grid = Grid_density

av_inc_z_dens = []

inc_grid_proj = []
z_grid_proj = []

inc_z_dens_Grid = np.zeros((len(min_z_grid),len(min_inc_grid)))

for z in range(len(min_z_grid)):
    for inc in range(len(min_inc_grid)):
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z])
        inc_grid_proj = np.append(inc_grid_proj,min_inc_grid[inc]+9)
        
        
        bin_data = []
        for d in range(len(dens_grid)):
            if Grid_z[d] == min_z_grid [z]:
                if Grid_inc[d] == min_inc_grid[inc]:
                    bin_data = np.append(bin_data,dens_grid[d])
        av_inc_z_dens = np.append(av_inc_z_dens,np.mean(bin_data))
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z])
        inc_grid_proj = np.append(inc_grid_proj,min_inc_grid[inc]+11)
        av_inc_z_dens = np.append(av_inc_z_dens,np.mean(bin_data))
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z]-0.25)     
        inc_grid_proj = np.append(inc_grid_proj,min_inc_grid[inc]+9)
        av_inc_z_dens = np.append(av_inc_z_dens,np.mean(bin_data))
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z]-0.5)
        inc_grid_proj = np.append(inc_grid_proj,min_inc_grid[inc]+9)
        av_inc_z_dens = np.append(av_inc_z_dens,np.mean(bin_data))
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z]-0.75)     
        inc_grid_proj = np.append(inc_grid_proj,min_inc_grid[inc]+9)
        av_inc_z_dens = np.append(av_inc_z_dens,np.mean(bin_data))
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z]-0.25)     
        inc_grid_proj = np.append(inc_grid_proj,min_inc_grid[inc]+11)
        av_inc_z_dens = np.append(av_inc_z_dens,np.mean(bin_data))
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z]-0.5)
        inc_grid_proj = np.append(inc_grid_proj,min_inc_grid[inc]+11)
        av_inc_z_dens = np.append(av_inc_z_dens,np.mean(bin_data))
        z_grid_proj = np.append(z_grid_proj,min_z_grid[z]-0.75)     
        inc_grid_proj = np.append(inc_grid_proj,min_inc_grid[inc]+11)
        av_inc_z_dens = np.append(av_inc_z_dens,np.mean(bin_data))

        inc_z_dens_Grid[z][inc] = np.mean(bin_data)

av_az_z_dens = []
az_grid_proj = []
z_grid_proj_az = []

az_z_dens_Grid = np.zeros((len(min_z_grid),len(min_az_grid)))

for z in range(len(min_z_grid)):
    for az in range(len(min_az_grid)):
        z_grid_proj_az = np.append(z_grid_proj_az,min_z_grid[z])
        az_grid_proj = np.append(az_grid_proj,min_az_grid[az])
        bin_data = []
        for d in range(len(dens_grid)):
            if Grid_z[d] == min_z_grid [z]:
                if Grid_az[d] == min_az_grid[az]:
                    bin_data = np.append(bin_data,dens_grid[d])
        av_az_z_dens = np.append(av_az_z_dens,np.mean(bin_data))

        z_grid_proj_az = np.append(z_grid_proj_az,min_z_grid[z]-0.25)     
        az_grid_proj = np.append(az_grid_proj,min_az_grid[az])
        av_az_z_dens = np.append(av_az_z_dens,np.mean(bin_data))
        z_grid_proj_az = np.append(z_grid_proj_az,min_z_grid[z]-0.5)
        az_grid_proj = np.append(az_grid_proj,min_az_grid[az])
        av_az_z_dens = np.append(av_az_z_dens,np.mean(bin_data))
        z_grid_proj_az = np.append(z_grid_proj_az,min_z_grid[z]-0.75)     
        az_grid_proj = np.append(az_grid_proj,min_az_grid[az])
        av_az_z_dens = np.append(av_az_z_dens,np.mean(bin_data))
        az_z_dens_Grid[z][az] = np.mean(bin_data)

##############################################################################################################################

################################  Plotting #################################################3

plt.scatter(z_grid_proj,inc_grid_proj,c=av_inc_z_dens,cmap=cmap,s=500,marker='s')
plt.xlim(-7,16)
plt.ylim(0,180)
plt.savefig('z_inc_proj.png')
plt.close()

plt.scatter(z_grid_proj_az,az_grid_proj,c=av_az_z_dens,cmap=cmap,s=500,marker='s')
plt.xlim(-7,16)
plt.ylim(0,360)
plt.savefig('z_az_proj.png')
plt.close()

