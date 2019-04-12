#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 09:00:29 2018

@author: lorenzo

Calculate average distance over all simulations for 1aps molecule CG potestio simulation at timestep = 0

This will be used as the benchmark target distance

"""
import os
import numpy as np

prot = "1aps"
os.chdir("/home/lorenzo.signorini/tirocinio/simulations/"+prot)

#######################
# loop through all runs
######################

l_of_lengths =[]

for runfolder in os.listdir():
    if runfolder.startswith("run"):
        print(runfolder)
        run_number = runfolder.lstrip("run_").rstrip("0").rstrip("_")
        os.chdir("/home/lorenzo.signorini/tirocinio/simulations/"+prot+"/"+runfolder+"/traj")
        traj=open("traj_"+run_number+".vtf")
        rl = traj.readlines()

        firstatom = False
        lastatom = False
        firststep = False
        n=0

        for line in rl:
            # extract coordinates of first and last atom when at the right line
            if firststep == True:
                n+=1
                if firstatom == True:
                    x = [float(i) for i in line.split()]
                    firstatom = False
                    
                if n == 98:
                    y = [float(i) for i in line.split()]                
              

            # check line of the timestep 0:
            if (line.startswith("t") and firststep == True):
                print("FINITO")
                firststep = False
                break
            
            if (line.startswith("t") and firststep == False): 
                print("CISIAMO")
                print(line)
                firststep = True
                firstatom = True
        
        # calculate distance
        dist= np.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2)
        l_of_lengths.append(dist)

avg = 0
for i in l_of_lengths:
    avg += i
avg = avg/len(l_of_lengths)

# print on file
os.chdir("/home/lorenzo/tirocinio/whole-atom-sMD")

f = open("avg_max_dist", "w")
header = "# OUTPUT OF mean_final_distance_calculator.py\n# contains list of all distances between first and last bead at timestep = 0 for 1aps for every cg simulation, as a benchmarck for the distance we would like to reach"
list_to_str = ",".join([str(i) for i in l_of_lengths])
f.write(header+"\n\n"+"average distance:    "+str(avg)+"\n\n"+"# Values:\n"+list_to_str)
f.close()



  
        
