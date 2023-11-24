# -*- coding: utf-8 -*-
"""
Differential Cross-section Analysis v2.0

Lukas Kostal, 18.11.2023, ICL
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.constants as sc
import scipy.optimize as so
import scipy.signal as ss


# gaussian function for fitting
def gauss(x, A, mu, sigma):
    y = A * np.exp(-(x - mu)**2 / (2 * sigma**2))    
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
def scat(tht_deg, Ei):
    tht = tht_deg / 180 * np.pi
    Ei *= sc.e * 1e3
    Ef = Ei / (1 + Ei / (sc.m_e * sc.c**2) * (1 - np.cos(tht)))
    Ef /= sc.e * 1e3
    return Ef


# function to calcaulte theoretical compton scattering
def cross(tht_deg, Ei):
    tht = tht_deg / 180 * np.pi
    Ei *= sc.e * 1e3
    Er = 1 / (1 + Ei / (sc.m_e * sc.c**2) * (1 - np.cos(tht))) 
    re = sc.physical_constants['classical electron radius'][0]
    cs = 1/2 * re**2 * Er**2 * (Er + 1/Er - np.sin(tht)**2)
    return cs


# function to load data
def load_mes(mes_num, mes_src, mes_tht, tg=True):
    # generate filepath for data with or without target
    if tg == True:
        filepath = f'Data/Mes{mes_num}_{mes_src}_{mes_tht:03}_tg.csv'
    if tg == False:
        filepath = f'Data/Mes{mes_num}_{mes_src}_{mes_tht:03}_ntg.csv'
    
    # load the data as a dataset
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
    E = E[~np.isnan(E)]
    N = N[~np.isnan(N)]
    
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
N_sim_min = 0

# factor of sigmas at which reduced chi sqaured is calculated
# take factor to get FWHM
f_sig = np.sqrt(2 * np.log(2))

# starting index of angles to include
ign = 2

# target distance from source in m
d_targ = 176 * 1e-3
derr_targ = 1 * 1e-3

# target mass in kg
m_targ = 84.31 * 1e-3
merr_targ = 0.01 *1e-3

# target molar mass in kg mol^-1 and proton no
Mm_targ = 26.981539 * 1e-3
Z_targ = 13

# detector aperture radius in m
r_det = 10 * 1e-3
rerr_det = 1 * 1e-3

# detector distance in m
d_det = 168 * 1e-3
derr_det = 1 * 1e-3

# load effective source activity in Bq
A_src, Aerr_src = np.loadtxt('Load/activity.csv', delimiter=',', unpack=True, skiprows=1)

# load the angles and associated errors in deg
geom = np.loadtxt('Load/geometry.csv', delimiter=',', unpack=True, skiprows=1)
tht_arr = (180 - geom[0, :]).astype(int)
tht_min = 180 - geom[2, :]
tht_max = 180 - geom[1, :]
tht_err = (tht_max - tht_min) / 2

# load reference peaks as dataframe
df = pd.read_csv('load/reference.csv', comment='#')

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

# arrays to hold peak yield and reduced chi sqaured from gaussian peak fit
yld_arr = np.zeros(len(tht_arr))
yld_err = np.zeros(len(tht_arr))
chi_arr = np.zeros(len(tht_arr))

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
    cmb = 1
    N_tg_cmb = N_tg
    h_tg_cmb = h_tg
    
    # whole max no of counts in bin is less than specified combine bins
    # note this is without subtracting the background
    while np.amax(N_tg_cmb) < N_mes_min:
        h_tg_cmb, N_tg_cmb = cmb_bin(h_tg, N_tg, cmb)
        cmb += 1
        
    h_tg_cmb, N_tg_cmb = cmb_bin(h_tg, N_tg, cmb)
    h_ntg_cmb, N_ntg_cmb = cmb_bin(h_ntg, N_ntg, cmb)
            
    # subtract the background rate and propagate error
    h = h_tg_cmb
    N = N_tg_cmb - N_ntg_cmb
    Nerr = np.sqrt(N_tg_cmb + N_ntg_cmb)
    
    # perform energy calibration
    E = poly(h, *pcal)
    
    # calculate expected energy of scattered gamma
    E0 = scat(tht, Ei)
    
    # find maximum no of counts in bin with at least 2 points in neighbourhood
    pks, _ = ss.find_peaks(N, width=2)
    N_max = np.amax(N[pks])
    
    # initial guess for preliminary curve fitting
    ig = np.array([N_max, E0, 10])

    # specify bounds on fitting parameters ±0.2 within the expected peak
    bnds_min = (N_max * 0.8, E0 * 0.8, 3)
    bnds_max = (N_max * 1.2, E0 * 1.2, 100)

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
    
    # calcualte yield of peak as full area of gaussian divided by measurement time
    yld = popt[0] * popt[2] * np.sqrt(2 * np.pi) / t_tg[-1]
    ylderr = yld * np.sqrt((perr[0] / popt[0])**2 + (perr[2] / popt[2])**2) / t_tg[-1]

    # write yield and reduced chi sqaured to arrays
    yld_arr[i] = yld
    yld_err[i] = ylderr
    chi_arr[i] = chi

    # array of energies for plotting
    E_plt = np.linspace(E[0], E[-1], 100)
    
    # parameters for plotting individual Compton peaks without saving
    plt.figure(1)
    plt.title('Compton Scattering Peak \n' \
              + f'measurement: Mes{mes_num}_{mes_src}_{tht:03}, ' \
              + f'$\chi^2=${chi:.2f}', pad=40)
        
    plt.xlabel(r'energy $E_\gamma$ (keV)')
    plt.ylabel(r'counts $N$ (unitless)')
    plt.rc('grid', linestyle=':', color='black', alpha=0.8)
    plt.grid()
    
    plt.errorbar(E, N, yerr=Nerr, fmt='.', c='blue', label='data')
    plt.plot(E_plt, gauss(E_plt, *popt), c='orangered', label='fit', zorder=np.inf)
    
    plt.legend(loc=(0.3, 1.02), ncol=2)
    plt.show()


# arrays to hold peak yield and reduced chi sqaured from gaussian peak fit
yld_sim_arr = np.zeros(len(tht_arr))
yld_sim_err = np.zeros(len(tht_arr))
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
    
    yld = np.sum(N_sim) / Ntot
    ylderr = np.sqrt(np.sum(Nerr_sim**2)) / Ntot
    
    # write yield and reduced chi sqaured to arrays
    yld_sim_arr[i] = yld
    yld_sim_err[i] = ylderr
    chi_sim_arr[i] = chi
    
    if tht < 20:
        yld_sim_arr[i] = np.nan
        yld_sim_err[i] = np.nan
        chi_sim_arr[i] = np.nan

# calcualte total no of e in target
Ne = m_targ / Mm_targ * sc.N_A * Z_targ
Ne_err = Ne * merr_targ / m_targ

# calcualte photon flux at target
flx = A_src / (4 * np.pi * d_targ**2)
flx_err = np.sqrt( (1 / (4 * np.pi * d_targ**2))**2 * Aerr_src**2 \
                + (3 * A_src / (4 * np.pi * d_targ**3))**2 * derr_targ**2)

# calculate detector solid angle
sa = np.pi * r_det**2 / d_det**2
sa_err = 2 * np.pi * np.sqrt( (r_det / d_det**2)**2 * rerr_det**2 \
                           + (r_det**2 / d_det**3)**2 * derr_det**2 )

# calcualte differnetial cross section from experimental yield
cs_arr = yld_arr / (Ne * flx * sa)
cs_err = cs_arr * np.sqrt( (yld_err / yld_arr)**2 + (Ne_err / Ne)**2 \
                         + (flx_err / flx)**2 + (sa_err / sa)**2)

# calcualte equivalent flux for simulation data
flx_sim = 1 / (4 * np.pi * d_targ**2)
flx_sim_err = np.sqrt((3 / (4 * np.pi * d_targ**3))**2 * derr_targ**2)

# calcualted differential cross section from simulated data
cs_sim_arr = yld_sim_arr / (Ne * flx_sim * sa)
cs_sim_err = cs_sim_arr * np.sqrt( (yld_sim_err / yld_sim_arr)**2 + (Ne_err / Ne)**2 \
                         + (flx_sim_err / flx_sim)**2 + (sa_err / sa)**2)

# note account for 
    
# calcualte theretical differential cross section
cs_exp_arr = cross(tht_arr, Ei)

# calcualte maximum and minimum theoretical cross section from systematic errors
cs_min = cross(tht_max, Ei - Eierr)
cs_max = cross(tht_min, Ei + Eierr)

# clacualte expected error on theoretical cross section
cs_exp_err = (cs_max - cs_min) / 2
    
# ignore the specified initial measurements at low angles
cs_arr = cs_arr[ign:]
cs_err = cs_err[ign:]
cs_sim_arr = cs_sim_arr[ign:] 
cs_sim_err = cs_sim_err[ign:] 
cs_exp_arr = cs_exp_arr[ign:]
cs_exp_err = cs_exp_err[ign:]
tht_arr = tht_arr[ign:]
tht_err = tht_err[ign:]


# combine errors from expected and observed variables for reduced chi sqaured fit
cs_err_chi = np.sqrt(cs_err**2 + cs_exp_err**2)

# calculate reduced chi sqaured value
chi = calc_chi(cs_arr, cs_err_chi, cs_exp_arr, 0)


# combine errors on theoretical and simulated cross section
cs_sim_err_chi = np.sqrt(cs_sim_err**2 + cs_exp_err**2)

# calcualte reduced chi sqaured value for simulated data
chi_sim = calc_chi(cs_sim_arr, cs_sim_err_chi, cs_exp_arr, 0)

# array of angles for plotting
tht_plt = np.linspace(0, 130, 100)

# parameters for plotting Compton scattering differential crossection
# note x axis units are bar 1b = 10^-28 m^2
plt.figure(1)
plt.title(f'Compton Scattering Differential Cross Section \n dataset: Mes{mes_num}, $\chi^2_\\nu=$ {chi:.2f}', pad=40)
plt.xlabel(r'scattering angle $\theta$ ($^{\circ}$)')
plt.ylabel(r'differential cross section $\frac{d \sigma}{d \Omega}$ ($b \; sr^{-1}$)')
plt.rc('grid', linestyle=':', color='black', alpha=0.8)
plt.grid()

plt.errorbar(tht_arr, cs_arr * 1e28, xerr=tht_err, yerr=cs_err * 1e28, fmt='.', c='blue', capsize=4, label='data')
plt.errorbar(tht_arr, cs_sim_arr * 1e28, xerr=tht_err, yerr=cs_sim_err * 1e28, fmt='.', c='green', capsize=4, label='simulation')
plt.plot(tht_plt, cross(tht_plt, Ei) * 1e28, c='orangered', label='theory', zorder=0)

plt.ylim(0, 0.08)

plt.legend(loc=(0.3, 1.02), ncol=2)

# save the plot
plt.savefig(f'Output/Crossection_mes{mes_num}.png', dpi=300, bbox_inches='tight')
plt.show() 

