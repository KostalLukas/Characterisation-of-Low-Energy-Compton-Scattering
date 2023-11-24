# -*- coding: utf-8 -*-
"""
Energy Calibration Analysis v3.1

Lukas Kostal, 14.11.2023, ICL
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as so
import scipy.signal as ss
import scipy.interpolate as si
import scipy.odr as sr
import itertools as it
import datetime


# gaussian function for fitting
def gauss(x, A, mu, sigma):
    y = A * np.exp(-((x - mu) / sigma)**2 / 2)    
    return y


# gaussian function for fitting with ODR
def odr_gauss(B, x):
    A, mu, sigma = B
    y = A * np.exp(-((x - mu) / sigma)**2 / 2)
    return y


# polynomial function for fitting
def poly(x, *coef):
    if type(x) == np.float64:
        y = 0
    else:
        y = np.zeros(len(x))
    for i in range(0, len(coef)):
        y += coef[i] * x**i
    
    return y


# polynomial function for fitting with odr
def odr_poly(coef, x):
    y = np.zeros(len(x))
    for i in range(0, len(coef)):
        y += coef[i] * x**i
    return y


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
    chi = np.sum(chi_arr) / (len(chi_arr) - dof)
    
    return chi


# dataset number to use for calibration
mes_num = 3

# minimum no of events for combining bins
N_min = 1000

# order of calibration polynomial
cal_ord = 3

# factor of sigmas at which reduced chi sqaured is calculated
f_sig = np.sqrt(2 * np.log(2))

# load reference peaks as dataframe
df = pd.read_csv('Load/reference.csv', comment='#')

# slice to array of sources and energies
src = df.iloc[:, 0].to_numpy()
E = df.iloc[:, 1].to_numpy()
Eerr = df.iloc[:, 2].to_numpy()

# sort the sources with the energies
i_srt = np.argsort(src)
src = src[i_srt]
E = E[i_srt]
Eerr = Eerr[i_srt]

# get array of unique reference sources so remove duplicates
src_ref, i_rpt, n_rpt = np.unique(src, return_index=True, return_counts=True)

# get the number of reference sources and reference peaks
n_ref = len(src_ref)
nn_ref = np.sum(n_rpt)

# array to append arrays of energies for reference peaks
E_ref = np.zeros(n_ref, dtype=object)
Eerr_ref = np.zeros(n_ref, dtype=object)

# loop over all reference sources write the peak energies to array
for i in range(0, n_ref):
    E_ref[i] = E[i_rpt[i]:i_rpt[i] + n_rpt[i]]
    Eerr_ref[i] = Eerr[i_rpt[i]:i_rpt[i] + n_rpt[i]]
    
# load the background measurement
t_bg, U_bg, n_bg, N_bg, g_bg = load_mes(mes_num, 'Bg')

# calculate scaled channel no by accounting for gain 
h_bg = n_bg / -g_bg

# calcualte total background measurement time and rate with error
T_bg = np.amax(t_bg)
R_bg = N_bg / T_bg
Rerr_bg = np.sqrt(N_bg) / T_bg

# array of reference peak energies and error in keV
E_pks = np.zeros(nn_ref)
Eerr_pks = np.zeros(nn_ref)

# array of scaled channel no from peak fit and error
h_pks = np.zeros(nn_ref)
herr_pks = np.zeros(nn_ref)

# array of reduced chi squared values for fitted peaks
chi_pks = np.zeros(nn_ref)

# loop over all reference sources with s as counting parameter
s=0
for i in range(0, n_ref):
    
    # load the data
    t, U, n, N, g = load_mes(mes_num, src_ref[i])
    
    # calculate scaled channel no by accounting for gain
    h = n / -g
    
    # calculate total measurement time and scale background and error
    T = np.amax(t)
    N_bg = R_bg * T
    Nerr_bg = Rerr_bg * T
    
    # prepare for combining bins
    cmb = 1
    N_cmb = N
    h_cmb = h
    
    # whole max no of counts in bin is less than specified combine bins
    while np.amax(N_cmb) < N_min:
        h_cmb, N_cmb = cmb_bin(h, N, cmb)
        cmb += 1
    h_cmb, N_cmb = cmb_bin(h, N, cmb)
    h_bg_cmb, N_bg_cmb = cmb_bin(h_bg, N_bg, cmb)

    # interpolate background with cubic spline for resampling
    bg_cs = si.CubicSpline(h_bg_cmb, N_bg_cmb)
    
    # resample background from cubic spline
    N_bg_rsmp = bg_cs(h_cmb)
    
    # subtract the background rate and propagate error
    h = h_cmb
    N = N_cmb - N_bg_rsmp
    Nerr = np.sqrt(N_cmb + np.abs(N_bg_rsmp))
    
    # minimum scaled channel below which peaks are not detected
    h_min = 0.08 * np.amax(h)
    i_min = np.argmin(np.abs(h - h_min))
    
    # find peaks and get their indices
    pks, _ = ss.find_peaks(N[i_min:], width=4)
    pks += i_min
    
    # no of reference peaks for the source
    n_pks = len(E_ref[i])
    
    # if no of found peaks is less than no of reference peaks raise exception
    if len(pks) < n_pks:   
        raise Exception('Error: insufficient number of peaks found')
    
    # sort the peaks by height and take the required no of peaks starting from heighest
    pks = pks[np.argsort(N[pks])[-n_pks:]]
    
    # manual correction for 2nd Na22 peak 
    if src_ref[i] == 'Na22':
        pks[0] = np.argmin(np.abs(h - 620))
    
    # array of optimised parameters, error on parameters, slice indices for each peak
    popt_arr = np.zeros(n_pks, dtype=object)
    perr_arr = np.zeros(n_pks, dtype=object)
    idx_arr = np.zeros(n_pks, dtype=object)
    
    # loop over each peak to fit
    for j in range(0, n_pks):
        
        # specify bounds on fitting parameters ±0.2 within the detected peak
        bnds_min = (N[pks[j]] * 0.8, h[pks[j]] * 0.8, 0)
        bnds_max = (N[pks[j]] * 1.2, h[pks[j]] * 1.2, 100)

        # initial parameteers from detecting peaks        
        ig = np.array([N[pks[j]], h[pks[j]], 10])
        
        # perform preliiminary gaussian fit over entire spectrum
        popt, pcov = so.curve_fit(gauss, h, N, p0=ig, bounds=(bnds_min, bnds_max))
        
        # find indices at specified factor of sigma below and above peak
        idx_min = np.argmin(np.abs(h - (popt[1] - f_sig * popt[2])))
        idx_max = np.argmin(np.abs(h - (popt[1] + f_sig * popt[2]))) +1
        
        # slice spectrum to given factor of sigma
        h_fit = h[idx_min:idx_max]
        N_fit = N[idx_min:idx_max]
        Nerr_fit = Nerr[idx_min:idx_max]
        
        # perform final gaussian fit over part of spectrum
        popt, pcov = so.curve_fit(gauss, h_fit, N_fit, p0=popt, sigma=Nerr_fit, absolute_sigma=True)
        
        # estimate error on fitting parameters from covariance matrix
        perr = np.sqrt(np.diag(pcov))
        
        # write potinised parameters, error on parameters and slicing indices to arrays
        popt_arr[j] = popt
        perr_arr[j] = perr
        idx_arr[j] = np.array([idx_min, idx_max])
        
        # calculate expected no of counts from gaussian fit
        N_exp = gauss(h_fit, *popt)
        
        # calcualte reduced chi sqaured value
        chi = calc_chi(N_fit, Nerr_fit, N_exp, 3)
    
        # write reference peak energy, scaled channel no from fit, reduced chi sqaured
        E_pks[s] = E_ref[i][j]
        h_pks[s] = popt[1]
        Eerr_pks[s] = Eerr_ref[i][j]
        herr_pks[s] = perr[1]
        chi_pks[s] = chi
        
        # increment counting parameter
        s +=1
        
    # title for plotting    
    tit = f'{src_ref[i]}, $g=${g:.2f}, $cmb=${cmb}, $\chi^2_\\nu=${chi_pks[s-1]:.1f}'
    
    # loop over each peak to append reduced chi squared value to title
    for j in range(1, n_pks):
        tit += f', {chi_pks[s-j-1]:.1f}'
    
    # plot each curve fit without saving
    plt.title(tit, pad=40)
    plt.xlabel(r'scaled MCA channel $h = \frac{n}{-g}$')
    plt.ylabel(r'counts $N$ (unitless)')
    plt.rc('grid', linestyle=':', color='black', alpha=0.8)
    plt.grid()
    
    plt.errorbar(h, N, yerr=Nerr, c='blue', label='data')
    plt.plot(h[pks], N[pks], '.', ms=8, c='red', zorder=8, label='peak')
    plt.axvline(x=h_min, ls='--', c='red', zorder=9, label='cutoff')

    # loop over each peak to plot gaussian fit
    for j in range(0, n_pks):
        h_plt = h[idx_arr[j][0]:idx_arr[j][1]]
        if j == 0:
            plt.plot(h_plt, gauss(h_plt, *popt_arr[j]), c='orangered', zorder=10, label='fit')
        else:
            plt.plot(h_plt, gauss(h_plt, *popt_arr[j]), c='orangered', zorder=10)
            
    plt.legend(loc=(0.1, 1.05), ncol=4)
        
    plt.show()

# array of 0 initial guess coefficients to specify polynomial order
ig = np.zeros(cal_ord+1)

# perform perliminary fit to get initial parameters for ODR
popt, pcov = so.curve_fit(poly, h_pks, E_pks, p0=ig)

# perform orthogonal distance regression with error in both variables
pmod = sr.Model(odr_poly)
pdat = sr.RealData(h_pks, E_pks, sx=herr_pks, sy=Eerr_pks)
podr = sr.ODR(pdat, pmod, beta0=popt)
pout = podr.run()

# get optimised parameters and error from orthogonal distance regression
popt = np.array(pout.beta)
perr = np.array(pout.sd_beta)

# array to hold expected scaled channel no from fit
h_exp = np.zeros(len(h_pks))

# loop over each reference peak
for i in range(0, len(h_pks)):
    
    # polynomial function with fit parameters offset by reference peak energy
    poly_inv = lambda x : poly(x, *popt) - E_pks[i]
    
    # use Newton-Raphson method to solve for expected scaled channel no
    h_exp[i] = so.newton(poly_inv, x0=h_pks[i])

# no of combinations of n object of 2 distinct types (2^n)
n_combs = 2**len(popt)

# create array of all possible combinations of +, -
combs = list(it.product([1, -1], repeat=len(popt)))
combs = np.array(combs)

# 2D array to hold expected scalled channel no for all possible combinations of fit parameters
h_exp_combs = np.zeros((len(h_pks), n_combs))

# loop over all possible combinations of errors on fit parameters
for i in range(0, n_combs):
    popt_combs = popt + combs[i] * perr

    # loop over each reference peak
    for j in range(0, len(h_pks)):
        
        # polynomial function with fit parameters offset by reference peak energy
        poly_inv = lambda x : poly(x, *popt_combs) - E_pks[j]
        
        # use Newton-Raphson method to solve for expected scaled channel no
        root, info = so.newton(poly_inv, x0=h_pks[j], disp=False, full_output=True)
        
        if info.converged and root > 0:
            h_exp_combs[j][i] = root
        else:
            h_exp_combs[j][i] = np.nan        

# take halve of maximum difference in combinations to be error in expected value
herr_exp = (np.nanmax(h_exp_combs, axis=1) - np.nanmin(h_exp_combs, axis=1)) / 2

# combine error in measured and expected scaled channel no in quadrature for chi calcualtion
herr_chi = np.sqrt(herr_pks**2 + herr_exp**2)

# calculate reduced chi sqaured for polynomial calibration
chi = calc_chi(h_pks, herr_chi, h_exp, len(popt))

# array of scaled channel no for plotting
h_plt = np.linspace(0, np.amax(h_pks), 100)

# 2D array to hold energies for all possible combinations of fit parameters
E_fit_combs = np.zeros((len(h_plt), n_combs))
E_res_combs = np.zeros((len(h_pks), n_combs))

# loop over all possible combinations of errors on fit parameters
# for each calcualte expected energy
for i in range(0, n_combs):
    popt_combs = popt + combs[i] * perr
    E_fit_combs[:, i] = poly(h_plt, *popt_combs)
    E_res_combs[:, i] = poly(h_pks, *popt_combs)

# take highest and lowest energies for plotting uncertainty region
E_fit_min = np.amin(E_fit_combs, axis=1)
E_fit_max = np.amax(E_fit_combs, axis=1)

# get residuals from polynomial calibration
E_res = E_pks - poly(h_pks, *popt)

# calculate error in residuals from error in fitting parameters
Eerr_res = np.sqrt(Eerr_pks**2 + (np.ptp(E_res_combs, axis=1) / 2)**2)

# get timestamp for the calibration
ts = datetime.datetime.now()
ts = ts.replace(microsecond=0)

# write calibration coefficients to .csv file
with open('Load/calibration.csv', 'w') as cal:
    print('order, coefficient, uncertainty', file=cal)
    
    for i in range(0, cal_ord+1):
        print(f'{i:.0f}, {popt[i]}, {perr[i]}', file=cal)
    
    print(f'# dataset: Mes{mes_num}', file=cal)
    print(f'# rcs value: {chi:2f}', file=cal)
    print(f'# timestamp: {ts}', file=cal)

# parameters for plotting calibration curve
plt.figure(1)
plt.title(f'Calibration Curve dataset: Mes{mes_num}, order: {cal_ord}, $\chi^2_\\nu=$ {chi:.2f}', pad=40)
plt.xlabel(r'scaled ADC channel $h = \frac{n}{-g}$')
plt.ylabel(r'energy $E_\gamma$ (keV)')
plt.rc('grid', linestyle=':', color='black', alpha=0.8)
plt.grid()

# plot the datapoints and uncertainty
plt.errorbar(h_pks, E_pks, xerr=herr_pks, yerr=Eerr_pks, fmt='.', c='blue', capsize=4, label='data')
plt.plot(h_plt, poly(h_plt, *popt), c='orangered', label='polynomial')
plt.fill_between(h_plt, E_fit_min, E_fit_max, color='royalblue', alpha=0.3, label='uncertainty')

plt.legend(loc=(0.1, 1.05), ncol=3)

# save the plot
plt.savefig(f'Output/Calibration_mes{mes_num}.png', dpi=300, bbox_inches='tight')
plt.show() 

# parameters for plotting residuals
plt.figure(2)
plt.title(f'Calibration Residuals dataset: Mes{mes_num}, order: {cal_ord}, $\chi^2_\\nu=$ {chi:.2f}')
plt.xlabel(r'scaled ADC channel $h = \frac{n}{-g}$')
plt.ylabel(r'energy $E_\gamma$ (keV)')
plt.rc('grid', linestyle=':', color='black', alpha=0.8)
plt.grid()

# plot the datapoints and uncertainty
#plt.plot(h_pks, E_res, 'x', c='blue')
plt.errorbar(h_pks, E_res, xerr=herr_pks, yerr=Eerr_res, fmt='.', c='blue', capsize=4)

# save the plot
plt.savefig(f'Output/Residuals_mes{mes_num}.png', dpi=300, bbox_inches='tight')
plt.show() 
