# -*- coding: utf-8 -*-
"""
Simulation Root Analysis v2.0

Lukas Kostal, 19.11.2023, ICL
"""


import numpy as np
import matplotlib.pyplot as plt
import uproot as ur
import glob as gb


# simulation number to convert
sim_num = 5

# no of bins for histogramming data 
n_bin = 1024

# total no of particles per run
Npart = 1

# total no of runs per dataset
Nrun = 2e8

# factor to account for restricted momentum vector
fmomRest = 0.5

# get a sorted list of all simulation datasets with given simulation number
sim = gb.glob(f'Simulation/output/Sim{sim_num}**.root')
sim.sort()

# get no of found simulation datasets
n_sim = len(sim)

# pring no of simulation files found
print(f'found {n_sim} root files from simulation: Sim{sim_num}')
print()

# loop over all simulation datasets
for i in range(0, n_sim):
    
    # get simulation dataset name
    sim_nam = sim[i].split('/')[-1][:-5]
    
    # print simulation root dataset being converted
    print(f'processing dataset: {sim_nam}.root')
    
    # load simulation root dataset
    ds_root = ur.open(sim[i])

    try: 
        # load events from root dataset    
        events = ds_root['events;1']
    except:
        # if get exception
        print(f'could not load events;1 in {sim_nam}.root')
    else:
        
        # load event no, total energy in keV and momentum vectors in keV
        Nevnt = events['Nevnt'].array()
        Etot = events['Etot'].array()
        Xmom = events['Xmom'].array()
        Ymom = events['Ymom'].array()
        Zmom = events['Zmom'].array()
    
        # find maxima for setting bin range
        Etot_max = np.amax(Etot)
        mom_max = np.amax(np.array([Xmom[:], Ymom[:], Zmom[:]]))
        
        # create bins for histogramming
        Etot_bin = np.linspace(0, Etot_max, n_bin)
        mom_bin = np.linspace(0, mom_max, n_bin)
        
        # arrays of bin centers
        Etot_cbin = Etot_bin[:-1] + np.diff(Etot_bin) / 2
        mom_cbin = mom_bin[:-1] + np.diff(mom_bin) / 2
        
        # histogramm simulated data
        Etot_hist, _ = np.histogram(Etot, Etot_bin)
        Xmom_hist, _ = np.histogram(Xmom, mom_bin)
        Ymom_hist, _ = np.histogram(Ymom, mom_bin)
        Zmom_hist, _ = np.histogram(Zmom, mom_bin)
    
        # column titles for csv dataset
        tits = 'event no, total energy (keV), E frequency, ' \
             + 'momentum component (keV), x frequency, y frequency, z frequency, ' \
             + 'particles per run, runs per simulation, restricted momentum direction'
    
        # ensure event number Nevnt is numpy array
        Nevnt = np.array(Nevnt)
        
        # get no of rows for each type of column
        n_Nevnt = len(Nevnt)
        n_Etot = len(Etot_cbin)
        n_mom = len(mom_cbin)
        
        # get max no of rows to be written in csv dataset
        n_max = np.amax(np.array([n_Nevnt, n_Etot, n_mom]))
        
        # write calibration coefficients to .csv file
        with open(f'Load/{sim_nam}.csv', 'w') as ds_csv:
            
            # print titles
            print(tits, file=ds_csv)
            
            # loop over max no of rows to be written
            for j in range(0, n_max):
                
                # string with event no data
                if j < n_Nevnt:
                    Nevnt_str = f'{Nevnt[j]}'
                else:
                    Nevnt_str = ''
                    
                # string with total energy data
                if j < n_Etot:
                    Etot_cbin_str = f'{Etot_cbin[j]}'
                    Etot_hist_str = f'{Etot_hist[j]:.0f}'
                else:
                    Etot_cbin_str = ''
                    Etot_hist_str = ''
                    
                # string with momentum components data
                if j < n_mom:
                    mom_cbin_str = f'{mom_cbin[j]}'
                    Xmom_hist_str = f'{Xmom_hist[j]:.0f}'
                    Ymom_hist_str = f'{Ymom_hist[j]:.0f}'
                    Zmom_hist_str = f'{Zmom_hist[j]:.0f}'
                else:
                    mom_cbin_str = ''
                    Xmom_hist_str = ''
                    Ymom_hist_str = ''
                    Zmom_hist_str = ''
                
                # string with simulation setup data
                if j == 0: 
                    Npart_str = f'{Npart:.0f}'
                    Nrun_str = f'{Nrun:.0f}'
                    momRest_str = f'{fmomRest}' 
                else:
                    Npart_str = ''
                    Nrun_str = ''
                    momRest_str = ''
                
                # combine strings to a row to be written to the csv file
                dats = Nevnt_str + ',' + Etot_cbin_str + ',' + Etot_hist_str + ',' \
                     + mom_cbin_str + ',' + Xmom_hist_str + ',' + Ymom_hist_str + ',' + Zmom_hist_str + ',' \
                     + Npart_str + ',' + Nrun_str + ',' + momRest_str
                
                # write the row to the csv file
                print(dats, file=ds_csv)
            
            # parameters for plotting simulated energy spectrum
            plt.figure(1)
            plt.title(f'Simulation Energy Sepctrum dataset:{sim_nam}.root')
                
            plt.xlabel(r'energy $E_\gamma$ (keV)')
            plt.ylabel(r'counts $N$ (unitless)')
            plt.rc('grid', linestyle=':', color='black', alpha=0.8)
            plt.grid()
            
            plt.plot(Etot_cbin, Etot_hist, c='blue')
            plt.ylim(0, np.amax(Etot_hist[:-1]) * 1.1)
            
            plt.show()

# print status update
print()
print('done')

