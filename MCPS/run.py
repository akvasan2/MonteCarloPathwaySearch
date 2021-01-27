#################### Configuration file for MCPS algorithm ############################
######################################################################################

import numpy as np
import math as m
from random import seed
from random import randint
import multiprocessing as mp
from multiprocessing import Pool
from joblib import Parallel, delayed
from datetime import datetime
import matplotlib.pyplot as plt
import MCPS as TS
import logging

################ Input parameters ################

data=np.loadtxt('Input_Files/input.dat') #### dataset taken from docking
output_file='Output_Files/transition_search.dat' ###### location for output paths from algorithm 

log_file='log.log' ####### log file to indicate errors and time to run algorithm

z_ind=0 # index for z
inc_ind=1 # index for inclination
az_ind=2 # index for azimuthal
energy_ind=3 # index for energy

top_z=np.max(data[:,0]) #topmost z to start search 
bott_z=np.min(data[:,0]) # bottommost z to end search
z_step=1 ### step size in z. 
inc_step=18 ### step size in inclination. 
az_step=18 ### step size in azimuthal. 

num=10 ### number of MCPS trajectories to obtain.  Check when convergence occurs to tune this parameter. 

proc = 2 ### number of processors to use (running is done in parallel)

##############################################create z bins to use in algorithm#########################################################

print("Creating z bins")

z_num=int((top_z-bott_z)/z_step)
N_zbins=z_num
z_list={}

for i in range(0,z_num):
    z_i=top_z-z_step*i
    z_list[i]=[]
    for f in range(len(data)):
        if (data[int(f)][z_ind]<(z_i) and data[int(f)][z_ind]>=(z_i-z_step)):
            z_list[i]=np.append(z_list[i],f)

#################################### run algorithm using specified number of processors (proc) ###############################

print("Starting MCPS")

log=open(log_file,'w')

try:
    d = datetime.now()
    starttime = datetime.now()
    log.write(f"Began process:{starttime} \n")
    parallel_searches=Parallel(n_jobs=proc)(delayed(TS.Transition_Search)(data,z_list,output_file,N_zbins,inc_step,az_step,z_ind,inc_ind,az_ind,energy_ind,itera) for itera in range(num-1))
    
    endtime = datetime.now()
    log.write(f"Finished process successfully {endtime}")

except:
    log.write(f"Job was terminated")
log.close()

