######################## Monte Carlo Based Pathway Search Algorithm ##############################

######################## Contributors: Archit Vasan and Nandan Haloi ############################
######################## Affiliation: University of Illinois at Urbana-Champaign########################
####################### Emails: akvasan2@illinois.edu, nhaloi2@illinois.edu, ######################

############### Code description: Method to use docking data to obtain a pathway which is connected (has no large jumps)#############################

################################# Connected in: z, inclination angle, and azimuthal angle space.  The values for these variables need to be previously calculated.
################################# The z-coordinates are used as a primary coordinate to define the pathway progressioni while the angles restrict exploration of the path.
################################# Moves are attempted using Metropolis criterion.

####################################################################################################################################################

import numpy as np
import math as m
from random import seed
from random import randint
import multiprocessing as mp
from multiprocessing import Pool
from joblib import Parallel, delayed
from datetime import datetime
import matplotlib.pyplot as plt
import csv


##### procedure used to standardize energies ######

def standardize_energy(energy_data):

	mean=np.mean(energy_data)
	std=np.std(energy_data)
	return [(e-mean)/float(std) for e in energy_data  ]

####### procedure to obtain allowable frames for the next iteration################

def allowed_frames(data,z_list,memory_array,prev_frame,inc_ind,az_ind,inc_step,az_step,z_i,indep=0):

	allowed_fr=[]
	allowed_fr_z_bin={}

##################### angle values from previous bin ###########################
	
	prev_inc=data[int(prev_frame)][inc_ind]
	prev_az=data[int(prev_frame)][az_ind]

###################### boundaries for angles in next step ############################ 

	inc_top=prev_inc+inc_step
	inc_bott=prev_inc-inc_step

######### need to take care of periodicity effects for azimuthal angle

	az_top=prev_az+az_step
	az_top_per=az_top - 360 ###### if az_top is above 360, this variable shifts the top boundary to the periodic image

	az_bott=prev_az-az_step
	az_bott_per=360+az_bott ###### if az_bot is below 0, this variable shifts the bottom boundary to the periodic image

################# Create the list of allowable frames using boundaries on angular space #######################

############## Create list using frames in the same z bin 

	for f in z_list[z_i]:
		n_f=int(f)
		next_inc=data[n_f][inc_ind]
		next_az=data[n_f][az_ind]


		########## Move in multiple angles in one step. only check frame if its memory_array value is 0 (not previously accepted) ####################
		if (indep==0 and memory_array[n_f]==0):

                        ######### if azimuthal bottom boundary is less than 0:
			if ( az_bott < 0 and next_inc>inc_bott and next_inc<inc_top and\
				( ( next_az>az_bott_per and next_az<360) \
				or (next_az<az_top and next_az > 0) ) ):

				allowed_fr=np.append(allowed_fr,n_f)
				allowed_fr_z_bin[n_f]=z_i

                        ####### if azimuthal top boundary is greater than 360 
			elif (az_top > 360 and next_inc>inc_bott and next_inc<inc_top and \
				( (next_az < az_top_per and next_az > 0) \
				or ( next_az > az_bott and next_az < 360 )) ):

				allowed_fr=np.append(allowed_fr,n_f)
				allowed_fr_z_bin[n_f]=z_i
			
			###### if az top boundary is < 360 and bott boundary >0, don't worry about periodicity
			elif (az_top < 360 and az_bott > 0 and next_inc>inc_bott and \
				next_inc<inc_top and \
				next_az>az_bott and next_az<az_top):

				allowed_fr=np.append(allowed_fr,n_f)
				allowed_fr_z_bin[n_f]=z_i

############ compute similar search for next z bin.  

	for f in z_list[int(z_i+1)]:
		n_f=int(f)
		next_inc=data[n_f][inc_ind]
		next_az=data[n_f][az_ind]

		if(memory_array[n_f]==0):
			if ( az_bott < 0 and next_inc>inc_bott and next_inc<inc_top and\
				( ( next_az>az_bott_per and next_az<360) \
				or (next_az<az_top and next_az > 0) ) ):
				allowed_fr=np.append(allowed_fr,n_f)
				allowed_fr_z_bin[n_f]=z_i+1

			elif (az_top > 360 and next_inc>inc_bott and next_inc<inc_top and \
				( (next_az < az_top_per and next_az > 0) \
				or ( next_az > az_bott and next_az < 360 )) ):
				allowed_fr=np.append(allowed_fr,n_f)
				allowed_fr_z_bin[n_f]=z_i+1

			elif (az_bott < 360 and az_bott > 0 and next_inc>inc_bott and \
				next_inc<inc_top and \
				next_az>az_bott and next_az<az_top):

				allowed_fr=np.append(allowed_fr,n_f)
				allowed_fr_z_bin[n_f]=z_i+1

	return allowed_fr,allowed_fr_z_bin ########### return both the frames as well as the z bin each frame belongs to

############### Procedure to perform actual Monte Carlo Based Pathway Search #############################

def Transition_Search(data,z_list,output_file,N_zbins,inc_step,az_step,z_ind,inc_ind,az_ind,energy_ind,iteration):

############### Define variables  and initialize arrays used in rest of code ##############################

	N_frames=len(data) ### Total poses in docking data
	memory_array=np.zeros(N_frames) #### array used to check if frame has already been accepted
	total_en_unstand=data[:,energy_ind]    
	en_stand=standardize_energy(total_en_unstand)   #### Standardize all total energies to use in algorithm

	z_i = 0 ############ start in 0th (top) z bin
	acc_frames=[]  ######## initialize array holding accepted frames
	acc_energies=[] ######## initialize array holding enery of acc frame

############ Pick random initial frame in topmost z bin ###########################

	seed(iteration)  ####### seed using user provided iteration
	rand=randint(0,len(z_list[0])-1)   ###### pick a random position in the beginning to start

	start_fr=int(z_list[0][rand])
	start_en=en_stand[start_fr]

	acc_frames=np.append(acc_frames,start_fr)  ####### add random frame to accepted frames array

############# While loop which performs pathway search until the last primary coordinate (z) bin is reached ###################

	step_num=0  ###### step along path identified so far

	while (z_i<N_zbins-1):

		prev_frame=acc_frames[int(step_num)] ########## previous frame identified in path
		prev_en=en_stand[int(prev_frame)]

################ Use above allowed_frames procedure to create allowed frames for next step of algorithm ##################

		allowed_frames_proc=allowed_frames(data,z_list,memory_array,prev_frame,inc_ind,az_ind,inc_step,az_step,z_i,indep=0) 
		allowed_fr=allowed_frames_proc[0]
		allowed_fr_z_bin=allowed_frames_proc[1]
		allowed_fr_num=len(allowed_fr)

		####### if no frames are available, termainate pathway search ########

		if (allowed_fr==[]):
			return 0

		accepted=0
		rejected=0
		bad_energy=0

################ Perform Metropolis search until a pose is accepted #############################

		while (accepted==0):
			######## if no frames allowed to be checked anymore, terminate pathway #######
			if(rejected>=allowed_fr_num or allowed_fr_num==[]):
				return 0
			####### choose random frame from allowed_fr (fr_try) to check acceptance 

			index=randint(0,int(len(allowed_fr)-1))
			fr_try=int(allowed_fr[index])
			en_try=en_stand[fr_try]

			met_crit=m.exp(prev_en-en_try) #Metropolis relation between previous energy and new one.

			###### automatically accept if met_crit >=1
			#####  if accepted, add to memory array so it wont be checked again, update z_i to new z bin, increase step_num
			if met_crit>=1:
				accepted+=1
				acc_frames=np.append(acc_frames,fr_try)
				z_i=allowed_fr_z_bin[fr_try]
				memory_array[fr_try]=1
				step_num+=1

			######### otherwise spit random number between 0, 1. if metrop> random number, accept

			else:
				rand=np.random.uniform(0,1)
				if met_crit>rand:
					accepted+=1
					acc_frames=np.append(acc_frames,fr_try)
					z_i=allowed_fr_z_bin[fr_try]
					memory_array[fr_try]=1
					step_num+=1
				else:
					it=np.where(allowed_fr==fr_try)
					allowed_fr=np.delete(allowed_fr,it) ##### remove from list of allowed poses
					memory_array[fr_try]=1 #### removing from memory if rejected as well right now
					rejected+=1
					continue
############### once a path is accepted write to output file ########################
	f=open(output_file,'a+',newline='')
	writer=csv.writer(f,delimiter=' ',quotechar='|', quoting=csv.QUOTE_MINIMAL)
	writer.writerow(acc_frames)
	return np.array(acc_frames)
