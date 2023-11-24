# -*- coding: utf-8 -*-
"""
Source Activity Analysis v2.0

Lukas Kostal, 18.11.2023, ICL
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.constants as sc
import scipy.optimize as so


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


# function to calcaulte theore
def cross(tht_deg, Ei):
    tht = tht_deg / 180 * np.pi
    Ei *= sc.e * 1e3
    Er = 1 / (1 + Ei / (sc.m_e * sc.c**2) * (1 - np.cos(tht))) 
    re = sc.physical_constants['classical electron radius'][0]
    cs = 1/2 * re**2 * Er**2 * (Er + 1/Er - np.sin(tht)**2)
    return cs
    

# function to load data
def load_mes(mes_num, src):
    # load the data as a dataset
    df = pd.read_csv(f'Data/Mes{mes_num}_{src}.csv')
    
    # slice into arrays
    t = df.iloc[:, 0].to_numpy()
    U = df.iloc[:, 4].to_numpy()
    N = df.iloc[:, 1].to_numpy()
    n = df.iloc[:, 2].to_numpy()
    g = df.iloc[0, -1]
    
    # remove nans from shorter columns
    t = t[~np.isnan(t)]
    U = U[~np.isnan(U)]
    n = n[~np.isnan(n)]
    N = N[~np.isnan(N)]
    
    return t, U, n, N, g


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
mes_num = 6

# source used in measurement
mes_src = 'Cs137'

# no of bins to combine
cmb = 0

# detector aperture radius in m
r_det = 10 * 1e-3
rerr_det = 1 * 1e-3

# detector distance in m
d_det = 194 * 1e-3
derr_det = 1 * 1e-3

# specified activity of the source in Bq and time ago when measured in days
A_spec = 3.7 * 1e6
Aerr_spec = 0.1 * 1e6 
t_spec = 1107
t_h = 11027

# factor of sigmas at which reduced chi sqaured is calculated
# take factor to get FWHM
f_sig = np.sqrt(2 * np.log(2))

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

# load target (tg) and no target (ntg) measurements
t, U, n, N, g = load_mes(mes_num, mes_src)
t_bg, U_bg, n_bg, N_bg, g_bg = load_mes(mes_num, "Bg")

# check if gain for tg and ntg measurements is same
if g != g_bg:
    raise Exception('Error: different source and background measurement gain')

# check if measurement time for tg and ntg measurement is same
if np.abs(t[-1] - t_bg[-1]) > 1:
    raise Exception('Error: different source and background measurement time')

# calculate scaled channel numbers
h = n / -g
h_bg = n_bg / -g_bg

if cmb > 0: 
    h, N = cmb_bin(h, N, cmb)
    h_bg, N_bg = cmb_bin(h_bg, N_bg, cmb)
        
# subtract the background rate and propagate error
N -= N_bg

Nerr = np.sqrt(N + N_bg)

# perform energy calibration
E = poly(h, *pcal)

# initial guess for preliminary curve fitting
ig = np.array([np.amax(N), Ei, 10])

# specify bounds on fitting parameters ±0.2 within the expected peak
bnds_min = (np.amax(N) * 0.8, Ei * 0.8, 3)
bnds_max = (np.amax(N) * 1.2, Ei * 1.2, 100)

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
yld = popt[0] * popt[2] * np.sqrt(2 * np.pi) / t[-1]
yld_err = yld * np.sqrt((perr[0] / popt[0])**2 + (perr[2] / popt[2])**2)

# instead just take activity to be total no of counts
yld = np.sum(N) / t[-1]
yld_err = np.sqrt(np.sum( Nerr**2 )) / t[-1]
 
# calculate effective activity of source
# note this includes detector peak efficiency
A = 4 * yld * d_det**2 / r_det**2
Aerr = 4 * np.sqrt( (d_det / r_det)**4 * yld_err**2 \
                 + (2 * yld * d_det / r_det**2)**2 * derr_det**2 \
                 + (2 * yld * d_det**2 / r_det**3)**2 * rerr_det**2 )

# calculate expected activity of the source
A_exp = A_spec * np.exp(-t_spec/t_h)    
Aerr_exp = A_spec * Aerr_spec / A_exp

# caluclate approximate detector efficiency and error
eff = A / A_exp
eff_err = eff * np.sqrt((Aerr / A)**2 + (Aerr_exp / A_exp)**2)

# print the effective activity into .csv file
with open('Load/activity.csv', 'w') as act_csv:
    print('effective activity (Bq), error in effective activity (Bq)', file=act_csv)
    print(f'{A}, {Aerr}', file=act_csv)

# print the calculated activity and detector efficiency
print(f'effective activity  = {A} ± {Aerr}')
print(f'detector efficiency = {eff:.4g} ± {eff_err:.4g}')

# array of energies for plotting
E_plt = np.linspace(E[0], E[-1], 100)

# parameters for plotting source activity peak
plt.figure(1)
plt.title('Source Activity Peak \n' \
          + f'measurement: Mes{mes_num}_{mes_src}, ' \
          + f'$\chi_\\nu^2=${chi:.2f}', pad=40)
    
plt.xlabel(r'energy $E_\gamma$ (keV)')
plt.ylabel(r'counts $N$ (unitless)')
plt.rc('grid', linestyle=':', color='black', alpha=0.8)
plt.grid()

plt.errorbar(E, N, yerr=Nerr, fmt='.', c='blue', label='data')
plt.plot(E_plt, gauss(E_plt, *popt), c='darkorange', label='fit', zorder=np.inf)

plt.legend(loc=(0.3, 1.02), ncol=2)

# save the plot
plt.savefig(f'Output/Activity_mes{mes_num}.png', dpi=300, bbox_inches='tight')
plt.show() 

