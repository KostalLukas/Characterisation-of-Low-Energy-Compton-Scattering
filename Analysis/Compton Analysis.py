# -*- coding: utf-8 -*-
"""
Compton Scattering Analysis v2.0

Lukas Kostal, 14.11.2023, ICL
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.constants as sc
import scipy.optimize as so


# gaussian function for fitting
def gauss(x, A, mu, sigma):
    y = A * np.exp(-((x - mu) / sigma)**2 / 2)    
    return y


# gaussian function for fitting with ODR
def odr_gauss(B, x):
    A, mu, sigma = B
    y = A * np.exp(-((x - mu) / sigma)**2 / 2)
    return y


# polynomial function for calibration
def poly(x, *coef):
    if type(x) == np.float64:
        y = 0
    else:
        y = np.zeros(len(x))
    for i in range(0, len(coef)):
        y += coef[i] * x**i
    
    return y


# function to calcaulte energy of scattered photons
# scattering angle tht in deg, scattered photon energy in keV
def scat(tht_deg, Ei):
    tht = tht_deg / 180 * np.pi
    Ei *= sc.e * 1e3
    Ef = Ei / (1 + Ei / (sc.m_e * sc.c**2) * (1 - np.cos(tht)))
    Ef /= sc.e * 1e3
    return Ef


# function to load measurement data
def load_mes(mes_num, mes_src, mes_tht, tg=True):
    # generate filepath for data with or without target
    if tg == True:
        filepath = f'Data/Mes{mes_num}_{mes_src}_{mes_tht:03}_tg.csv'
    if tg == False:
        filepath = f'Data/Mes{mes_num}_{mes_src}_{mes_tht:03}_ntg.csv'
    
    # load the data as a dataframe
    df = pd.read_csv(filepath, comment='#')
    
    # slice into arrays
    t = df.iloc[:, 0].to_numpy()
    U = df.iloc[:, 4].to_numpy()
    N = df.iloc[:, 1].to_numpy()
    n = df.iloc[:, 2].to_numpy()
    g = df.iloc[0, -1]
    
    # remove nans from shorter columns
    # also ignore last channel since it will contain overflow
    t = t[~np.isnan(t)]
    U = U[~np.isnan(U)]
    n = n[~np.isnan(n)][:-1]
    N = N[~np.isnan(N)][:-1]
    
    return t, U, n, N, g


# function to load simulation data
def load_sim(sim_num, sim_tht):
    # generate simulated dataset filepath
    filepath = f'Load/Sim{sim_num}_{sim_tht:03}.csv'
    
    # load the data as a dataframe
    df = pd.read_csv(filepath, comment='#')
    
    # slice into arrays
    Nevnt = df.iloc[:, 0].to_numpy()
    E = df.iloc[:, 1].to_numpy()
    N = df.iloc[:, 2].to_numpy()
    # note here can also load momentum components
    Npart = int(df.iloc[0, 7])
    Nrun = int(df.iloc[0, 8])
    fmomRest = float(df.iloc[0, 9])
    
    # remove nans from shorter columns
    # also ignore last channel since will contain overflow
    Nevnt = Nevnt[~np.isnan(Nevnt)]
    E = E[~np.isnan(E)][:-1]
    N = N[~np.isnan(N)][:-1]
    
    # calculate total no of simulated particles
    # correct for restricted momentum
    Ntot = Npart * Nrun / fmomRest
    
    return Nevnt, E, N, Ntot


# function to combine bins
def cmb_bin(n_arr, N_arr, cmb):
    cmb = int(cmb)
    n_arr = n_arr[(len(n_arr) % cmb):]
    N_arr = N_arr[(len(N_arr) % cmb):]
    
    N_arr = np.sum(N_arr.reshape(-1, cmb), axis=1)
    n_arr = np.mean(n_arr.reshape(-1, cmb), axis=1)

    return n_arr, N_arr


# function to get reduced chi squared
def calc_chi(o, sig, e, dof):
    with np.errstate(all='ignore'):
        chi_arr = (o - e)**2 / (sig)**2
        
    chi_arr = chi_arr[~np.isinf(chi_arr)]
    chi_arr = chi_arr[~np.isnan(chi_arr)]
    chi = np.sum(chi_arr) / (len(chi_arr) - dof)
    return chi

# measurement number to analyse
mes_num = 5

# source used in measurement
mes_src = 'Cs137'

# simulation number to analyse
sim_num = 5

# minimum no of counts for combining bins
N_mes_min = 200
N_sim_min = 12

# factor of sigmas at which reduced chi sqaured is calculated
f_sig = np.sqrt(2 * np.log(2))

# load the angles and associated errors in deg
geom = np.loadtxt('Load/geometry.csv', delimiter=',', unpack=True, skiprows=1)
tht_arr = (180 - geom[0, :]).astype(int)
tht_min = 180 - geom[2, :]
tht_max = 180 - geom[1, :]
tht_err = (tht_max - tht_min) / 2

# load reference peaks as dataframe
df = pd.read_csv('Load/reference.csv', comment='#')

# slice to array of sources and energies
src = df.iloc[:, 0].to_numpy()
E = df.iloc[:, 1].to_numpy()
Eerr = df.iloc[:, 2].to_numpy()

# find initial energy of source used and associated error
Ei = E[src == mes_src][0]
Eierr = Eerr[src == mes_src][0]

# load calibration file
cal = np.loadtxt('Load/calibration.csv', delimiter=',', unpack=True, skiprows=1)

# get calibration parameters and error
pcal = cal[1]

# arrays to hold final energy of gamma from fiting and reduced chi sqaured
Ef_arr = np.zeros(len(tht_arr))
Ef_err = np.zeros(len(tht_arr))
chi_arr = np.zeros(len(tht_arr))

# loop over all measured datasets
for i in range(0, len(tht_arr)):
    
    # angle of current measurement in deg
    tht = tht_arr[i]

    # load target (tg) and no target (ntg) measurements
    t_tg, U_tg, n_tg, N_tg, g_tg = load_mes(mes_num, mes_src, tht, tg=True)
    t_ntg, U_ntg, n_ntg, N_ntg, g_ntg = load_mes(mes_num, mes_src, tht, tg=False)

    # check if gain for tg and ntg measurements is same
    if g_ntg != g_tg:
        raise Exception('Error: different tg and ntg measurement gain')

    # check if measurement time for tg and ntg measurement is same
    if np.abs(t_tg[-1] - t_ntg[-1]) > 1:
        raise Exception('Error: different tg and ntg measurement time')

    # calculate scaled channel numbers
    h_tg = n_tg / -g_tg
    h_ntg = n_ntg / -g_ntg

    # prepare for combining bins
    cmb_mes = 1
    N_tg_cmb = N_tg
    h_tg_cmb = h_tg
    
    # if whole max no of counts in bin is less than specified combine bins
    # note this is without subtracting the background
    while np.amax(N_tg_cmb) < N_mes_min:
        h_tg_cmb, N_tg_cmb = cmb_bin(h_tg, N_tg, cmb_mes)
        cmb_mes += 1
    h_tg_cmb, N_tg_cmb = cmb_bin(h_tg, N_tg, cmb_mes)
    h_ntg_cmb, N_ntg_cmb = cmb_bin(h_ntg, N_ntg, cmb_mes)
            
    # subtract the background rate and propagate error
    h = h_tg_cmb
    N = N_tg_cmb - N_ntg_cmb
    Nerr = np.sqrt(N_tg_cmb + N_ntg_cmb)
    
    # perform energy calibration
    E = poly(h, *pcal)

    # calculate expected energy of scattered gamma
    E0 = scat(tht, Ei)
    
    # initial guess for preliminary curve fitting
    ig = np.array([np.amax(N), E0, 10])

    # specify bounds on fitting parameters ±0.2 within the expected peak
    bnds_min = (np.amax(N) * 0.8, E0 * 0.8, 3)
    bnds_max = (np.amax(N) * 1.2, E0 * 1.2, 100)

    # perfom preliminary gaussian fit
    popt, pcov = so.curve_fit(gauss, E, N, p0=ig, bounds=(bnds_min, bnds_max))
    
    # find indices at specified factor of sigma below and above peak
    idx_min = np.argmin(np.abs(E - (popt[1] - f_sig * popt[2])))
    idx_max = np.argmin(np.abs(E - (popt[1] + f_sig * popt[2]))) +1
    
    # slice spectrum to given factor of sigma
    E_fit = E[idx_min:idx_max]
    N_fit = N[idx_min:idx_max]
    Nerr_fit = Nerr[idx_min:idx_max]
    
    # perform final gaussian fit over part of spectrum
    popt, pcov = so.curve_fit(gauss, E_fit, N_fit, p0=ig, sigma=Nerr_fit, absolute_sigma=True)
    
    # estimate error on fitting parameters from covariance matrix
    perr = np.sqrt(np.diag(pcov))
    
    # expected no of counts from gaussian fit
    N_exp = gauss(E_fit, *popt)

    # calcualte reduced chi sqaured for the fit
    chi = calc_chi(N_fit, Nerr_fit, N_exp, 3)

    # write results to the arrays
    Ef_arr[i] = popt[1]
    Ef_err[i] = perr[1]
    chi_arr[i] = chi
    
    # array of energies for plotting
    E_plt = np.linspace(E_fit[0], E_fit[-1], 100)
    
    # parameters for plotting individual Compton peaks without saving
    plt.figure(1)
    plt.title('Compton Scattering Peak \n' \
              + f'measurement: Mes{mes_num}_{mes_src}_{tht:03}, ' \
              + f'$\chi^2=${chi:.2f}', pad=40)
        
    plt.xlabel(r'energy $E_\gamma$ (keV)')
    plt.ylabel(r'counts $N$ (unitless)')
    plt.rc('grid', linestyle=':', color='black', alpha=0.8)
    plt.grid()
    
    plt.errorbar(E, N, yerr=Nerr, fmt='.', c='blue', label='simulation')
    plt.plot(E_plt, gauss(E_plt, *popt), c='orangered', label='fit', zorder=np.inf)
    
    plt.legend(loc=(0.3, 1.02), ncol=2)
    plt.show()

# arrays to hold final energy of gamma from fiting and reduced chi sqaured
Ef_sim_arr = np.zeros(len(tht_arr))
Ef_sim_err = np.zeros(len(tht_arr))
chi_sim_arr = np.zeros(len(tht_arr))

# loop over all simualted datasets
for i in range(0, len(tht_arr)):
    
    # angle of current measurement in deg
    tht = tht_arr[i]

    # now also load simulated data
    Nevnt, E_sim, N_sim, Ntot = load_sim(sim_num, tht)
    
    # prepare for combining bins for simulated data
    cmb_sim = 1
    N_sim_cmb = N_sim
    E_sim_cmb = E_sim
    
    # if whole max no of counts in bin is less than specified combine bins
    while np.amax(N_sim_cmb) < N_sim_min:
        E_sim_cmb, N_sim_cmb = cmb_bin(E_sim, N_sim, cmb_sim)
        cmb_sim += 1
    E_sim, N_sim = cmb_bin(E_sim, N_sim, cmb_sim)

    # get error on simulated data
    Nerr_sim = np.sqrt(N_sim)
    
    # if error is 0 take it to be 0.1 instead
    Nerr_sim[Nerr_sim < 1] = 0.1 
    
    # fit at all angles expect 0
    if tht > 10:
        # calculate expected energy of scattered gamma
        E0 = scat(tht, Ei)
        
        # initial guess for preliminary curve fitting
        ig = np.array([np.amax(N_sim), E0, 10])
    
        # specify bounds on fitting parameters ±0.2 within the expected peak
        bnds_min = (np.amax(N_sim) * 0.8, E0 * 0.8, 5)
        bnds_max = (np.amax(N_sim) * 1.2, E0 * 1.2, 30)
    
        # perfom preliminary gaussian fit
        popt, pcov = so.curve_fit(gauss, E_sim, N_sim, p0=ig, bounds=(bnds_min, bnds_max))
        
        # find indices at specified factor of sigma below and above peak
        idx_min = np.argmin(np.abs(E_sim - (popt[1] - f_sig * popt[2])))
        idx_max = np.argmin(np.abs(E_sim - (popt[1] + f_sig * popt[2]))) +1
        
        # ensure enough datapoints for fitting
        if (idx_max - idx_min) < 4:
            idx_min -= 2
            idx_max += 2
        
        # slice spectrum to given factor of sigma
        E_fit = E_sim[idx_min:idx_max]
        N_fit = N_sim[idx_min:idx_max]
        Nerr_fit = Nerr_sim[idx_min:idx_max]
        
        # perform final gaussian fit over part of spectrum
        popt, pcov = so.curve_fit(gauss, E_fit, N_fit, p0=popt, sigma=Nerr_fit, absolute_sigma=True, maxfev=2000)
        
        # estimate error on fitting parameters from covariance matrix
        perr = np.sqrt(np.diag(pcov))
        
        # expected no of counts from gaussian fit
        N_exp = gauss(E_fit, *popt)

        # calcualte reduced chi sqaured for the fit
        chi = calc_chi(N_fit, Nerr_fit, N_exp, 3)
        
        # write fitted energy to arrays
        Ef_sim_arr[i] = popt[1]
        Ef_sim_err[i] = perr[1]
        chi_sim_arr[i] = chi
        
        # array of energies for plotting
        E_plt = np.linspace(E_fit[0], E_fit[-1], 100)
        
        # parameters for plotting individual Compton peaks without saving
        plt.figure(1)
        plt.title('Compton Scattering Peak \n' \
                  + f'measurement: Mes{mes_num}_{mes_src}_{tht:03}, ' \
                  + f'$\chi^2=${chi:.2f}', pad=40)
            
        plt.xlabel(r'energy $E_\gamma$ (keV)')
        plt.ylabel(r'counts $N$ (unitless)')
        plt.rc('grid', linestyle=':', color='black', alpha=0.8)
        plt.grid()
        
        plt.errorbar(E_sim, N_sim, yerr=Nerr_sim, fmt='.', c='green', label='data')
        plt.plot(E_plt, gauss(E_plt, *popt), c='orangered', label='fit', zorder=np.inf)
        
        plt.legend(loc=(0.3, 1.02), ncol=2)
        plt.show()
        
    # if theta is 0 ignore the measurement
    if tht < 20:
        Ef_sim_arr[i] = np.nan
        Ef_sim_err[i] = np.nan
        chi_sim_arr[i] = np.nan
    

# get array of theoretical energies of gamma after scattering
Et_arr = scat(tht_arr, Ei)

# for negative min angle take 0 since scattering is symmetric so would give 0 error in expected
tht_min[tht_min < 0] = 0

# get lowest and highest expected energy
Et_min = scat(tht_max, Ei - Eierr)
Et_max = scat(tht_min, Ei + Eierr)

# take halve difference between lowest and highest to be error on expected energy
Et_err = (Et_max - Et_min) / 2

# combine error in measurement and expected value for reduced chi sqaured
# note ignoring error in reference Ei since negligable
Ef_err_chi = np.sqrt(Ef_err**2 + Et_err**2)
Ef_sim_err_chi = np.sqrt(Ef_sim_err**2 + Et_err**2)

# calcualte reduced chi sqaured for Compton curve
chi = calc_chi(Ef_arr, Ef_err_chi, Et_arr, 0)
chi_sim = calc_chi(Ef_arr, Ef_sim_err_chi, Et_arr, 0)

# array of values of theta for plotting
tht_plt = np.linspace(np.amin(tht_arr), np.amax(tht_arr), 100)

# parameters for plotting Compton scattering curve
plt.figure(1)
plt.title(f'Compton Scattering Curve dataset: Mes{mes_num}, $\chi^2_\\nu=$ {chi:.2f}', pad=30)
plt.xlabel(r'scattering angle $\theta$ ($^{\circ}$)')
plt.ylabel(r'energy $E_\gamma$ (keV)')
plt.rc('grid', linestyle=':', color='black', alpha=0.8)
plt.grid()

# plot the datapoints and uncertainty
plt.errorbar(tht_arr, Ef_arr,xerr=tht_err, yerr=Ef_err, fmt='.', c='blue', capsize=4, label='data')
plt.errorbar(tht_arr, Ef_sim_arr,xerr=tht_err, yerr=Ef_sim_err, fmt='.', c='green', capsize=4, label='simulation')
plt.plot(tht_plt, scat(tht_plt, Ei), c='orangered', label='theory', zorder=0)

plt.plot(tht_arr, Ef_sim_arr, '.', c='green')

plt.legend(loc=(0.2, 1.02), ncol=3)

# save the plot
plt.savefig(f'Output/Compton_mes{mes_num}.png', dpi=300, bbox_inches='tight')
plt.show() 
