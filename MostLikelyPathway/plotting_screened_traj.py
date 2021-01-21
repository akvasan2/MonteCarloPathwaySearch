import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits import mplot3d
import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)
#######################################

########## Load in all screened trajectories, cluster1 trajectories, and cluster2 trajectories ########

t_file=open('Traj_Group_Data/screened_trajectories.dat')
transitions_tot=[]
with t_file as my_file:
    for line in my_file:
        myarray=np.fromstring(line, dtype=float, sep=' ')
        transitions_tot.append(myarray)
transitions=transitions_tot

t_file=open('Traj_Group_Data/cluster1_paths.dat')
transitions_tot=[]
with t_file as my_file:
    for line in my_file:
        myarray=np.fromstring(line, dtype=float, sep=' ')
        transitions_tot.append(myarray)
transitions_c1=transitions_tot

t_file=open('Traj_Group_Data/cluster2_paths.dat')
transitions_tot=[]
with t_file as my_file:
    for line in my_file:
        myarray=np.fromstring(line, dtype=float, sep=' ')
        transitions_tot.append(myarray)
transitions_c2=transitions_tot

######### Load in input data #################

data = np.loadtxt('../MCPS/Input_Files/input.dat')
#######################################


######################FOr each group of trajectories, make lists of the azimuthal, inclination angles#############################
######################################################################################################################################

transitions_inc = {}
transitions_az = {}
transitions_z = {}
for i in range(len(transitions_c2)):
    transitions_inc[i] = []
    transitions_az[i] = []
    transitions_z[i] = []
    for j in range(len(transitions_c2[i])):
        f = int(transitions_c2[i][j])
        transitions_z[i] = np.append(transitions_z[i], data[f][0])
        transitions_inc[i] = np.append(transitions_inc[i], data[f][1])
        transitions_az[i] = np.append(transitions_az[i], data[f][2])

transitions_inc_c1 = {}
transitions_az_c1 = {}
transitions_z_c1 = {}
for i in range(len(transitions_c1)):
    transitions_inc_c1[i] = []
    transitions_az_c1[i] = []
    transitions_z_c1[i] = []
    for j in range(len(transitions_c1[i])):
        f = int(transitions_c1[i][j])
        transitions_z_c1[i] = np.append(transitions_z_c1[i], data[f][0])
        transitions_inc_c1[i] = np.append(transitions_inc_c1[i], data[f][1])
        transitions_az_c1[i] = np.append(transitions_az_c1[i], data[f][2])

transitions_inc_c2 = {}
transitions_az_c2 = {}
transitions_z_c2 = {}
for i in range(len(transitions_c2)):
    transitions_inc_c2[i] = []
    transitions_az_c2[i] = []
    transitions_z_c2[i] = []
    for j in range(len(transitions_c2[i])):
        f = int(transitions_c2[i][j])
        transitions_z_c2[i] = np.append(transitions_z_c2[i], data[f][0])
        transitions_inc_c2[i] = np.append(transitions_inc_c2[i], data[f][1])
        transitions_az_c2[i] = np.append(transitions_az_c2[i], data[f][2])


########### Plot inclination vs z for each cluster #############

for i in range(len(transitions_c1)):
    plt.plot(transitions_z_c1[i],transitions_inc_c1[i],color='black')


for i in range(len(transitions_c2)):
    plt.plot(transitions_z_c2[i],transitions_inc_c2[i],color='forestgreen')

plt.ylim(0,180)
plt.savefig('Images/clusters_z_inc.png',bbox_inches='tight')
plt.close()

########### Plot azimuthal vs z for each cluster ################

transitions_c1_az_total =  []
transitions_c1_inc_total = []
transitions_c1_z_total = []

transitions_c2_az_total =  []
transitions_c2_inc_total = []
transitions_c2_z_total = []

for i in range(len(transitions_c1)):
    for j in range(len(transitions_c1[i])):
        if transitions_az_c1[i][j]>290:
            transitions_az_c1[i][j]-=360 ##### This is to deal with periodic trajectories
        transitions_c1_az_total = np.append(transitions_c1_az_total,transitions_az_c1[i][j])
        transitions_c1_z_total = np.append(transitions_c1_z_total,transitions_z_c1[i][j])
    plt.plot(transitions_z_c1[i],transitions_az_c1[i],color='black')

for i in range(len(transitions_c2)):
    for j in range(len(transitions_c2[i])): 
        if transitions_az_c2[i][j]<75:
            transitions_az_c2[i][j]+=360 ##### This is to deal with periodic trajectories
        transitions_c2_az_total = np.append(transitions_c2_az_total,transitions_az_c2[i][j])
        transitions_c2_z_total = np.append(transitions_c2_z_total,transitions_z_c2[i][j])

    plt.plot(transitions_z_c2[i],transitions_az_c2[i],color='forestgreen')

plt.ylim(0,360)
plt.savefig('Images/clusters_z_az.png',bbox_inches='tight')
plt.close()

