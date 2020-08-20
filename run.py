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

data=np.loadtxt('Input_Files/input.dat')
output_file='Output_Files/transition_search_ang18_solv_en.dat'
log_file='log.log'

top_z=np.max(data[:,0])
bott_z=np.min(data[:,0])
z_step=1
inc_step=18
az_step=18
small_dev_inc=18
small_dev_az=18
z_ind=0
inc_ind=1
az_ind=2
energy_ind=3

num=10000
proc = 16

z_num=int((top_z-bott_z)/z_step)
N_zbins=z_num
z_list={}

for i in range(0,z_num):
    z_i=top_z-z_step*i
    z_list[i]=[]
    for f in range(len(data)):
        if (data[int(f)][z_ind]<(z_i) and data[int(f)][z_ind]>=(z_i-z_step)):
            z_list[i]=np.append(z_list[i],f)


log=open(log_file,'w')

try:
    d = datetime.now()
    starttime = datetime.now()
    log.write(f"Began process:{starttime} \n")
    print("Test")
    parallel_searches=Parallel(n_jobs=proc)(delayed(TS.Transition_Search)(data,z_list,output_file,N_zbins,inc_step,az_step,small_dev_inc,small_dev_az,z_ind,inc_ind,az_ind,energy_ind,itera) for itera in range(num-1))
    
    endtime = datetime.now()
    log.write(f"Finished process successfully {endtime}")

except:
    log.write(f"Job was terminated")
log.close()

