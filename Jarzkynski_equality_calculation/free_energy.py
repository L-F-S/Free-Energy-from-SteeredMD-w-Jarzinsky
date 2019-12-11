#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 11:25:34 2018

@author: Lorenzo Federico Signorini


###############################################################################
##         ALGORITHM FOR CALCULATING THE FREE ENERGY                         ##
           DIFFERENCE FROM A FOLDED TO AN UNFOLDED
           STATE OF A PROTEIN USING JARZYNSKY EQUALITY
           ON STEERED MOLECULAR DYNAMICS SIMULATIONS
 
 For reference:
    . Jarzynski, Christopher. "Nonequilibrium equality for free energy differences."
    Physical Review Letters 78.14 (1997): 2690.
    . Sanghyun Park, Fatemeh Khalili-Araghi, Emad Tajkhorshid, and Klaus Schulten.
      "Free energy calculation from steered molecular dynamics simulations using jarzynski’s
    equality." The Journal of chemical physics, 2003

Requires as input, several steered molecular dynamics (SMD) simulations.
calculates the work done for every simulation,
then calcualtes free energy from the work, according to Jarzynsky's equality.

Works with NAMD simulatiosn file outputs.

------------------------------------------------------------------------------
V 5.3 (tipo) : simulazioni fatte in modo diverso  calcolo del lavoro rifatto,
ciclando su distanze incrementali.

update: trapezoidal rule for integral calculation

update: z-test for cheching normality of distirbution.

update: aggiunto il grafico delle stianze su tempo real vs expected (preso da v4)

update: Levato W_with_t, adesso

update: levata la correzione per smd delle vecchie simulazioni

update: cambio delle paths per ottenere le posizioni dell'smd atom, xke h
il nuovo tipo di simulazione.
update: cambiato il time, nel nuovo formato, ogni line ha il tempo in ps nella
colon

uptdate: calcolodel lavoro per ogni angstrom di distanza usando una matrice

todo: (per la 5.3) crtl+f todo . correggere il lavoro dividendo tutto per 100,

###############################################################################
#                                                                             #
#                                  MAIN                                       #
###############################################################################
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as ss

###############################################################################
# Define useful global variables

T=310 # in kelvin
kb = 0.0019872041 # kcal/mol*K
beta = 1/(T*kb) # in mol/Kcal
k = 7          # force constant.  kcal/(mol*A^2)

# utility varables for calculations
speed = "0005"  
actual_v= "0."+speed
v_angstrom_fsec = float(actual_v)/2 
n_simulations = 2 # not actually nsimulations, ma il numero che metot nel ciclo: 1001, 200, 200, 12

# Useful folders:
pdb_dir = "/home/lorenzo/home/lorenzo.signorini/tirocinio/inputs/"
smd_dir = "/home/lorenzo.signorini/tirocinio/out/"+speed
ana_dir = "/home/lorenzo/tirocinio/whole-atom-sMD/analyses/nuove_nuove_simu/"+speed

# set maximum distance at which to calculate work.
max_dist = 100 # sarebbe tipo 350 - 26 che è tipo una cosa sensata) xke la linea rdistance è 372.4
dist_0 =  30 # you know this for a FACT. BECAUSE YOU ENFORCED THIS IN THE SIMULATION.

def main(load_works = False, load_deltaF = False):
    print("Hello! wellcome to this calculation of the difference in\n\
          free energy from a folded state (at an end-to-end distance of "\
           , dist_0,"to an unfolded state (at an end-to-end distance of ", max_dist,\
           "from steered molecular dynamic simulation data of protein pulling\
           at speed: ", speed,"\n+++++++++++++++++++++++++++++++++++++++++++++++++")
   
    #0. Plot expected vs real distance vs time.
   
    plot_all_dist() 
   
    #1. Calculate Work:

    if load_works == False:
        W_with_d = calculate_all_works_with_matrix() 
        save_work(W_with_d)          
    elif load_works == True:
        W_with_d = upload_work()
        
    W_with_d = W_with_d.dropna()  # Nan from rows of experimenti che non ci sono
    print(len(W_with_d.index ), "experiments")
    print( W_with_d.shape[1], "position steps.")
    total_work = W_with_d[max_dist].dropna()
    print(total_work)
    
    # 1.5 Work Plots
    
    plot_all_w_wrt_dist(W_with_d)
    histo(total_work)
    # is_normal(total_work) normality check
    
    #2 calc free energy
    
    cumulant = use_cumulant(total_work)

    if load_deltaF == False:
        DeltaF_of_dist = calc_free_energy(W_with_d, cumulant)
         save_DeltaF(DeltaF_of_dist, cumulant) todo: deltaf è una lista, non ha .to_csv
    else:
        DeltaF_of_dist = upload_DeltaF(cumulant)
    
    print("FINAL ENERGY:", DeltaF_of_dist[-1])

    # 2.5 DeltaF PLOTS
    
    DeltaF_plot(DeltaF_of_dist, sorted(W_with_d.columns), cumulant)
    return #W_with_d,# DeltaF_of_dist

"""
###############################################################################
#                                                                             #
#                        1. CALCULATE WORK                                    #
###############################################################################

Calculates work applied on the system from t=o to t=finaltime,
by discretizing the integral of the original work

w(0->finale)= k*v* sum_(i=0->i=t_finale){[(x_(i+1)-xi)-v*t_i]* valore della tstep}  
(unità di misura: kcal/mol)

where xi is the distance between C_alpha N terminal and C_alpha C terminal
at time i  

todo: add description of trapezoid rulke discretization.
"""    

def extract_quantity_inside_sommatoria(x_ser, line, dist_i_1, trapezio = False):
    time = float(line.split()[0])*1000 # *1000 per rimetterlo in femtosecondi, perche' e' scritto in pico       
    deltat = 100 # grandezza timestep in fs
    i = time/float(deltat)
    # extract current position of nterm. aggiornato a v4
    x_tyr_i=[float(i) for i in line.split()[1:4]]        
    # calculate distance from nterm to cterm
    dist_i = np.sqrt((x_tyr_i[0]-x_ser[0])**2+(x_tyr_i[1]-x_ser[1])**2+(x_tyr_i[2]-x_ser[2])**2)    
    # calculate the quantity inside the sommatoria
    if trapezio == True:
        quantita = dist_i + dist_i_1 - 2*dist_0 - v_angstrom_fsec*(2*i-1)
        return quantita, dist_i, time     
    quantita = dist_i - dist_0 - v_angstrom_fsec*time
    print("step", i)
    print("quantits", quantita)
    print("chiamata a ", i)
    return quantita, dist_i, time

def extract_coordinates_of_atom_of_simulation(i, bound):
    if bound == True:
    # extract coordinates of Ca of N-term SER:
        pdbfile =open(pdb_dir+"aps_start_"+str(i)+".pdb")
        for line in pdbfile.readlines():
            if line.split()[1] == "5":    # l'indice del CA della Ser
                    x = [float(i) for i in line.split()[6:9]]
                    break
        pdbfile.close()
    
    else:
    # extract coordinates of bound atom:
        pdbfile =open(pdb_dir+"aps_start_"+str(i)+".pdb")
        for line in pdbfile.readlines():
            if line.split()[1] == "1530":    # l'indice del CA della Ser
                    x = [float(i) for i in line.split()[6:9]]
                    break
        pdbfile.close()        
       
    return x

def get_initial_dist(i):
    x_ser = extract_coordinates_of_atom_of_simulation(i, True)
    x_tyr_0 = extract_coordinates_of_atom_of_simulation(i, False)
    dist_0_eff = np.sqrt((x_tyr_0[0]-x_ser[0])**2+(x_tyr_0[1]-x_ser[1])**2+(x_tyr_0[2]-x_ser[2])**2)  #distance at first step. should be = dist_0
    return x_ser, dist_0_eff

def calculate_all_works_with_matrix():
    """ FOR EVERY SPEED, Creates a pd.Dataframe (or a np.array), that has
    as rows the single pulling experiments, and as columns the expected distance
    in $\AA$, from 30 to 352. The values inside the DataFrame are the calculated 
    actual Work done on the system, for that experiment, until that distance.
    Time does not come into consideration, since the actual time is always >= 
    t expected , for every distance d=d0+v*t(expected)"""
    
    os.chdir(smd_dir)
    
    "-1. Initialize matrix"
    Work_dataset = pd.DataFrame(columns = range(dist_0, max_dist+1), index = range(1, n_simulations), dtype = 'float64')
    print("Creating matrix of works done at speed: ", speed,".\nNumber of experiments: ", Work_dataset.shape[0],".\nNumber of datapoints:  ", Work_dataset.shape[1],"   min: ", dist_0, "max: ", max_dist)
    
    "0. initialize last_current_max_distance, to save things in"
    last_current_max_distance = dist_0
    
    for current_max_distance in range(dist_0, max_dist+1):
        "1. Calculate work until dist = current_max_distance"
        print("Calculating works for distance = ", current_max_distance, "for all simulations.\n")
        for j in range(1,n_simulations):
            "2. Iterate all simulations (rows) "
            print("trajectory", j)
            filename= "aps_"+speed+"_"+str(j)+"_smd"
            if os.path.exists(filename):           # check file actually existis
                
                
                """a. set initial contitions: calculate initial positions of
                ser_ca, tyr_ca. Calculate initial distance (should be = dist_0).
                Set last_current_max_distance = 0, for the first iteration.
                (this variable will be used to save computational time)"""
                x_ser, dist_i = get_initial_dist(j)
    
                """ b. Initialize cumsum, which will store the values inside the
                sum in the formula for the calculation of work"""
                cumsum = 0  # to calculate work
                f=open(filename)
                for line in f.readlines(): 
                    if not line.startswith("g1"):

                        print("cumsum",cumsum)
                        "c. update dist_i, and caculate quantita."
                        quantita, dist_i, time = extract_quantity_inside_sommatoria(x_ser, line, dist_i, trapezio = False)

                        if dist_i < last_current_max_distance:
                            "do not calculate work for these steps"
                            print('not calc')
                            continue
                        else:
                            print(dist_i)
                            if dist_i <= current_max_distance: 
                                cumsum += quantita  
                                print(cumsum)
                            else:
                                # nb sei sempre nel for sulle simulazioni
                                break 
                f.close()        

                
                """3. Insert Work done for experiment j, until 
                distance = current_max_distance"""
                
                if current_max_distance == dist_0:
                     Work_dataset[current_max_distance][j] = - cumsum * k * v_angstrom_fsec * 100  # col, row 100: n of fs per step.
                else:
                    Work_dataset[current_max_distance][j] = - ( - Work_dataset[current_max_distance-1][j]/( k *v_angstrom_fsec * 100) + cumsum) * k *v_angstrom_fsec * 100
        # else: non esiste: do nothing
        #   WARNING: THIS WILL YIELD NAN IN ROWS SENZA SIMULAZIONI. TO BE TAKEN INTO ACCOUNT NEI GRAFICI
        "4. update last_current_max distance, for next iteration"
        last_current_max_distance = current_max_distance
    return Work_dataset

"""
###############################################################################
#                                                                             #
#                       2. Free energy calculations                            #
###############################################################################
"""
def use_cumulant(total_work):
    mean =  np.mean(total_work.dropna())
    var = np.var(total_work.dropna())
    print("Total Work:\nMean:",mean, "Variance:", var)
    if beta*mean >= (beta**2/2)*var:
        print("Calculating exponential average using cumulant expansion")
        return True
    else:
        print("Calculating exponential average as is")
        return False
    
def calc_free_energy(Work_at_d, cumulant):
    print("Calculating FREE ENERGY")
    DeltaF_of_dist = []
    for dist in sorted(Work_at_d.columns):
        works = np.array(Work_at_d[dist])
        if cumulant == True:
            DeltaF = -(1/beta)*(-beta*(np.mean(works))+((beta**2)/2)*np.var(works))
#            third_term = - (beta**3/6)*(np.mean(last_work**3) - 3*np.mean(last_work**2)*np.mean(last_work) + 2*np.mean(last_work)**2)
        else:
            DeltaF = -(1/beta)*np.log(np.mean(np.exp(-beta*works)))
        print(DeltaF)    
        DeltaF_of_dist.append(DeltaF)
    return DeltaF_of_dist

"""
###############################################################################
#                                                                             #
#                             3. PLOTS                                        #
###############################################################################
"""

def is_normal(total_work):
    os.chdir(ana_dir)
    ss.probplot(total_work, plot = plt)
    plt.savefig("Q-Q_plot_total_work_"+speed+".png")
    return True

def plot_all_dist():
    
    print("Plotting distance/time graphs with expected vs observed trajectories,\n for all trajectories:")
    os.chdir(smd_dir)

    missing_indexes=[]
    for i in range(1,n_simulations):
        filename= "aps_"+speed+"_"+str(i)+"_smd"
        distances = []
        times = []
        if os.path.exists(filename): # check file exists (sometimes alcuni possono mancare xke è crashata la simulazione, o altro)
            x_ser, dist_j = get_initial_dist(i)
            f=open(filename)    
            for line in f.readlines(): 
                if not line.startswith("g1"):
                    quantita, dist_j, time = extract_quantity_inside_sommatoria(x_ser, line, dist_j)
                    if dist_j<= max_dist:
                        distances.append(dist_j)
                        times.append(time)
                    else:
                        break 
            f.close()
            plot_dist_wrt_time_of_traj(i, distances, times)
        else:
            missing_indexes.append(i)
    return

def histo(final_work):
    os.chdir(ana_dir)
    final_work = np.array(final_work.dropna())  # xke some rows are na 
    plt.figure(figsize=(30, 20))
    matplotlib.rc('xtick', labelsize=25) 
    matplotlib.rc('ytick', labelsize=25) 
    plt.title("Pulling work (Kcal/mol) from 30 to 372 $\AA$ at "+"0."+speed+"$\AA$/fs", fontsize = 35)
    plt.hist(final_work, bins = 100, color = "black")
    plt.savefig("Final_work_"+speed+".png")
    plt.close()
    return

def DeltaF_plot(DeltaF_of_dist, dist, cumulant):
    os.chdir(ana_dir)
    plt.figure(figsize=(30, 20))
    matplotlib.rc('xtick', labelsize=25) 
    matplotlib.rc('ytick', labelsize=25) 
    plt.plot(dist,DeltaF_of_dist, color='black', linewidth=0.9)
    plt.axhline(y=0, color='grey', linestyle='--', linewidth=0.9)
    plt.xlim(0)
    plt.xlabel("Distance ($\AA$)", fontsize = 30)
    plt.ylabel("DeltaF (Kcal/mol)" , fontsize = 30)
    if cumulant == True:
        cum = "cumulant_"
    else:
        cum = "exp_"
    plt.savefig("DeltaF_"+cum+speed+".png")
    plt.close()
    return

def  plot_dist_wrt_time_of_traj(i, distances, times):
    os.chdir(ana_dir+"/distance_plots")
    print("trajectory ", i)
    plt.figure()
    plt.plot(np.array(times)/1000, distances)
    lambdavt = [v_angstrom_fsec*t + dist_0 for t in times]
    plt.plot(np.array(times)/1000, lambdavt)
    plt.xlabel("time (ps)")
    plt.ylabel("distance ($\AA$)")
    plt.savefig("distance_wrt_time"+"_"+speed+"_traj"+str(i)+".png")
    plt.close()
    os.chdir(smd_dir)
    return

def plot_all_w_wrt_dist(W_with_d):
    os.chdir(ana_dir+"/w_vs_dist")
    print("Plotting Work as a function of distance, for every trajectory:")
    for i, row in W_with_d.iterrows():
        print("trajectory", i)
        plt.figure(figsize=(30, 20))    
        plt.plot(sorted(W_with_d.columns), W_with_d.loc[i])
        plt.savefig("W_wrt_dist_traj"+str(i)+".png")
        plt.close()
    return

"""
###############################################################################
#                                                                             #
#                              4. I/O                                         #
###############################################################################
"""

def save_DeltaF(DeltaF, cumulant):
    os.chdir(ana_dir)
    if cumulant == True:
        cum = 'cum_'
    else:
        cum = 'exp_'
    filename="DeltaF_"+cum+speed+".csv"
    DeltaF.to_csv(filename, sep = ',', header = True, index = True)
    return

def upload_DeltaF(cumulant):
    print("Uploading works from previous calculations")
    os.chdir(ana_dir)
    
    if cumulant == True:
        cum = 'cum_'
    else:
        cum = 'exp_'
    filename="DeltaF_"+cum+speed+".csv"
    DeltaF_of_dist = pd.read_csv(filename, index_col = 0)
    return DeltaF_of_dist

def save_work(W_with_d):
    os.chdir(ana_dir)
    filename="work_"+speed+".csv"
    W_with_d.to_csv(filename, sep = ',', header = True, index = True)
    return

def upload_work():
    print("Uploading works from previous calculations")
    os.chdir(ana_dir)
    filename="work_"+speed+".csv"
    W_with_d = pd.read_csv(filename, index_col = 0)
    W_with_d.columns = list(range(dist_0, max_dist+1))
    return W_with_d

# run program:
main()
