#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper: "Elucidating the Interplay Between Entropy-Driven and Patch-Mediated 
Bonding in Directing Nanoscale Assemblies"
Sample Code for Chain Size Calculation
Authors: Kireeti Akkunuri, Xiangyu Zhang, Thi Vo
Dated: September 2024

"""

import gsd.hoomd
import numpy as np
import freud.box
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
import os
import glob

# read trajectory from gsd file
pattern = 'sample_single_part*.gsd'
file_list = glob.glob(pattern)

# Sort the list of files in ascending order
file_list.sort()

# Create empty object to store data
file = []
positions_read = []
images_read = []
timesteps_read = []
boxes_read = []
charge_read=[]
body_read=[]

# initial cutoff frame range
# this is to avoid collecting data from overstretched bonds initially
# these values below have to correlate with the main run code
ndump = 5E3
init_cutoff = int((500*1E4)/(ndump))
dt = 2.5E-4
# init_cutoff = 1000 # frames



# Check if any files were found
if not file_list:
    print("No files found matching the pattern.")
else:
    # Loop through and open each file using gsd.hoomd.open
    for file_name in file_list:
        try:
            with gsd.hoomd.open(name=file_name, mode='rb') as f:
                lenf=len(f)
                print(lenf)
                for ind, frame in enumerate(f):
                    # if ind >= init_cutoff and ind%5==0:
                    time_tmp = frame.configuration.step/ndump
                    if time_tmp >= init_cutoff:
                        timesteps_read.append(time_tmp*ndump*dt)
                        positions_read.append(frame.particles.position)
                        images_read.append(frame.particles.image)
                        boxes_read.append(frame.configuration.box)
                    if ind==lenf-1:
                        charge_read=frame.particles.charge
                        body_read=frame.particles.body
                        
                    

        except Exception as e:
            print(f"Error opening {file_name}: {e}")
# del f, lenf
del file

####### Combine segments of trajectory #######################################

def remove_redundancies(timesteps_read, positions_read, images_read, 
                        boxes_read):
    # Create a dict to keep track of the last occurrence index of each timestep
    last_occurrence_indices = {}
    for index, timestep in enumerate(timesteps_read):
        last_occurrence_indices[timestep] = index

    # Extract the indices of the last occurrences
    retained_indices = sorted(last_occurrence_indices.values())

    # Retain only the items at the retained indices for all lists
    timesteps_read = [timesteps_read[i] for i in retained_indices]
    positions_read = [positions_read[i] for i in retained_indices]
    images_read = [images_read[i] for i in retained_indices]
    boxes_read = [boxes_read[i] for i in retained_indices]

    # Combine the lists into tuples and sort them by the timesteps
    combined = sorted(zip(timesteps_read, positions_read, images_read, boxes_read))
    
    # Unzip the combined list back into individual lists
    timesteps_read, positions_read, images_read, boxes_read = map(list, zip(*combined))

    return timesteps_read, positions_read, images_read, boxes_read

timesteps_read, positions_read, images_read, boxes_read = remove_redundancies(timesteps_read, positions_read, images_read, boxes_read)
    
last = len(timesteps_read)-1
Nsnaps = last+1

N_particles = len(positions_read[0])
N_monomers = N_particles/3
chains = charge_read
bodies = body_read
N_bodies = int(max(bodies)+1*3)
N_chains = int(max(chains)+1)
N_chains = int(N_monomers/(N_bodies/3))
# Insphere diameter of shape
sigma=40.0

del bodies

timesteps=[]
rg2_series = np.zeros([0,2])    # Time series of average rg2
rg2_all = np.zeros([0,N_chains])
SD_list = np.zeros([0,2])


t0 = int(timesteps_read[0]/ndump)

# Collect snapshot data at each frame and compute r_g2
for i in range(0,Nsnaps,1):

    box = freud.Box.from_box(boxes_read[i])
    positions_uw = box.unwrap(positions_read[i],images_read[i])

    j = 0       # chain index
    posn = 0    # along chain, of monomer
    chn_values, chn_counts = np.unique(chains,return_counts=True)
    rg2_list = np.zeros((0,1)) # rg2 of each body in a snapshot

    
    while j < N_chains:
        # Calculate only for monomers at indices 0, 3, 6, 9...
        var_tmp = np.var(positions_uw[posn:posn+chn_counts[j]:3,:],axis=0,dtype=np.float32)
        rg2_tmp = ((var_tmp[0]+var_tmp[1]+var_tmp[2]))
        rg2_list = np.append(rg2_list,rg2_tmp)
        posn = posn + chn_counts[j]    
        j += 1
        
    # Every single rg2 value
    rg2_all = np.append(rg2_all,[rg2_list],
                          axis=0)
    rg2_series=np.append(rg2_series,[[timesteps_read[i],np.mean(rg2_list)]],
                          axis=0)
    SD_list=np.append(SD_list,[[timesteps_read[i],np.std(rg2_list,ddof=1)/np.sqrt(N_chains)]],axis=0)

# Calculating averages over last 40% of trajectory    
avg_rg2 = np.mean(rg2_series[int(0.6*len(rg2_series)):len(rg2_series)-1,1])
avg_SD = np.std(rg2_series[int(0.6*len(rg2_series)):len(rg2_series)-1,1])
del rg2_all 
with open('r_g.txt','w') as f:
 	print(int(float(N_monomers)/float(N_chains)),round(avg_rg2,1),round(avg_SD,1),file=f)
