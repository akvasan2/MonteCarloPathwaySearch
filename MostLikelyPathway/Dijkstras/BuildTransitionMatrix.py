import numpy as np
import math as m
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

################################### Run transition matrix code #########################################################
####################D###################D##################D

def BuildTransitionMatrix(transitions, data, z_ind, inc_ind, az_ind, en_ind, beta, top_z, bott_z, z_step, inc_step, az_step, z_num, bins_z, bins_inc_ang, bins_az_ang, Met_mat_file, Neg_log_file):
    transitions_grids={}
    transitions_energies={}
    for i in range(len(transitions)):
            transitions_grids[i] = []
            transitions_energies[i] = []
            for j in range(len(transitions[i])):
                f = int(transitions[i][j])
                f_b4 = int(transitions[i][j-1])
                if data[f_b4][z_ind]<=bins_z[len(bins_z)-1] and data[f_b4][z_ind]>=bins_z[len(bins_z)-1]-z_step:
                    break
                it = 0
                #### Assign point along pathway to a particular grid 
                for z in range(len(bins_z)):
                    if data[f][z_ind]<=bins_z[z] and data[f][z_ind]>=bins_z[z]-z_step:
                        for inc in range(len(bins_inc_ang)):
                            if data[f][inc_ind]>=bins_inc_ang[inc] and data[f][inc_ind]<=bins_inc_ang[inc]+inc_step:
                                for az in range(len(bins_az_ang)):
                                    it = (z*len(bins_inc_ang)*len(bins_az_ang)) + inc*len(bins_az_ang) + az 
                                    if data[f][az_ind]>=bins_az_ang[az] and data[f][az_ind]<=bins_az_ang[az]+az_step:
                                        transitions_grids[i] = np.append( transitions_grids[i], it)
                                        transitions_energies[i] = np.append(transitions_energies[i],data[f][en_ind])
                                        break
    
    ##########################################################################################
    
    it=0
    sources=np.zeros(len(transitions))
    sinks=np.zeros(len(transitions))
    
    for i in range(len(transitions)):
        sources[i] = transitions_grids[i][0] 
        sinks[i] = transitions_grids[i][-1]
    
    np.savetxt('sources.dat',sources)
    np.savetxt('sinks.dat',sinks)
    #########################################################################################
    
    N_clusters=len(bins_z)*len(bins_inc_ang)*len(bins_az_ang)+1
    
    Met_matrix=np.zeros((N_clusters,N_clusters))
    
    for p in range(len(transitions_grids)):
        for f in range(len(transitions_grids[p])-1):
                c_i=int(transitions_grids[p][f])
                c_iplus1=int(transitions_grids[p][f+1])
                e_i=transitions_energies[p][f]#data[int(transitions[p][f])][en_ind]
                e_iplus1=transitions_energies[p][f+1]#data[int(transitions[p][f+1])][en_ind]
                Met_matrix[c_i][c_iplus1]+=np.exp((e_i-e_iplus1)/float(beta))
    
    np.savetxt(Met_mat_file,Met_matrix) ##### Outputting a matrix of each energy weighted transition
    print("Met matrix created")
    
    #########################################################################################
    P_matrix = np.zeros((N_clusters+2,N_clusters+2))
    neg_log_mat=np.zeros((N_clusters+2,N_clusters+2))####### we are adding a virtual node to the start and end of the matrix.  In order for this to be a squre matrix, it is added to rows+ columns
    
    for ci_0 in range(1,N_clusters):
        for cj_0 in range(1,N_clusters):
            ci=ci_0-1 # Indices associated with met matrix
            cj=cj_0-1 # Indices associated with met matrix
            if Met_matrix[ci][cj]!=0:
                P_matrix[ci_0][cj_0]=Met_matrix[ci][cj]
            if ci_0!=cj_0 and P_matrix[ci_0][cj_0]!=0:
                neg_log_mat[ci_0][cj_0]=-np.log(P_matrix[ci_0][cj_0])
    
    for ci_0 in range(1,N_clusters+1):
        ci=ci_0-1
        if ci in sources:
            neg_log_mat[0][ci_0]=-50
        if ci in sinks:
            neg_log_mat[ci_0][-1]=-50
    
    np.savetxt(Neg_log_file,neg_log_mat)
    print("Transition matrix created")
    #########################################################################################
    return neg_log_mat

