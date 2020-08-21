# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 17:25:58 2020

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.signal import argrelextrema

##### INPUT ######
data = np.genfromtxt(r'C:\Users\User\Documents\03_software\python\PyMS\testdata\peakFitter_testdata\data2.txt')
# Initial guesses for the parameters to fit:
# amplitudes, means and standard deviations plus a continuum offset.
guess_wo_A = [(4342.8, 0.5),
              (4363.8, 0.5),
              (4403.7, 0.5),
              (4426.55, 0.5),
              (4450.75, 0.5),
              (4465.6, 0.5),
              (4492.2, 0.5),
              (4505.3, 0.5)]

peak_type = 'lorentzian'

ub_width = 30
lb_width = 0

ub_position = 1
lb_position = -1

offset = 0
ub_offset = 10
lb_offset = 0

charge = 11

# m/z region for peak finding
region_of_interest = (4320,4520)

##### DONT CHANGE FROM HERE #####

def sum_duplicate_x(data):
    """
    Sum up all y values at the same x-value.
    Useful for curating data coming from MassLynx' copy spectrum list
    """
    x_old = data[0][0]
    y_old = data[0][1]
    
    x_new = []
    y_new = []
    
    for x, y in data[1:]:
        if x == x_old:
            y+= y_old
        else:
            x_new.append(x)
            y_new.append(y)
            y = 0
        x_old = x
        y_old = y
        
    return np.array([x_new,y_new]).T

def gaussian(x, A, x0, sig):
    return A*np.exp(-np.pi*(x-x0)**2/(sig**2))

def lorentzian(x, A, x0, sig):
    return A*sig**2/(sig**2 + (2*x-2*x0)**2)

def multi_peaks(x, *pars):

    offset = pars[-1]
    ps = offset
    # split parameter list into chunks of three elements (leaving out the offset)
    lst = [pars[i:i+3] for i in range(0,len(pars)-1,3)]
    for i in lst: 
        if peak_type == 'gaussian':
            ps += gaussian(x, *i)
        elif peak_type == 'lorentzian':
            ps += lorentzian(x, *i)
    return ps

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

data = sum_duplicate_x(data)
x, y = data.T
# index of values inside regin of itnerest
xinside_roi = np.where((x > region_of_interest[0]) & (x < region_of_interest[1]))

if guess_wo_A == None:
    x_maxima = x[argrelextrema(y, np.greater)]
    # index of values inside regin of itnerest
    inside_roi = np.where((x_maxima > region_of_interest[0]) & (x_maxima < region_of_interest[1]))
    x_maxima = x_maxima[inside_roi]
    guess_wo_A = list(zip(x_maxima, [0.5]*len(x_maxima)))

guess = []
#guess_wo_A = [guess_wo_A[i:i+2] for i in range(0,len(guess_wo_A),2)]
for g in guess_wo_A:
    A = y[find_nearest_idx(x, g[0])]
    guess.extend([A, *g])
guess.append(offset)

lb_mz = [x + lb_position for x in guess[1::3]]
ub_mz = [x + ub_position for x in guess[1::3]]

lb_width = [lb_width]*len(guess_wo_A)
ub_width = [ub_width]*len(guess_wo_A)

lb_amp = [0] * len(guess_wo_A)
ub_amp = [np.inf] * len(guess_wo_A)
# fit and calculate optimised parameters and covariance matrix
popt, pcov = optimize.curve_fit(multi_peaks,
                                x[xinside_roi], y[xinside_roi], guess,
                                bounds=(list(np.concatenate(list(zip(lb_amp, lb_mz, lb_width))).flat)+[lb_offset,],
                                       list(np.concatenate(list(zip(ub_amp, ub_mz, ub_width))).flat)+[ub_offset]),
                                method='trf')

# reshape the optimised gaussian parameters (w/o the offset)
popt_res = popt[:-1].reshape(int((len(popt)-1)/3),3).tolist()

# calc one standard deviation errors on the parameters 
perr = np.sqrt(np.diag(pcov))
# reshape the output
perr_res = perr[:-1].reshape(int((len(perr)-1)/3),3).tolist()

print('====PeakFitter====\n')
for idx, tup in enumerate(popt_res):
    print('Peak {}:'.format(idx))
    print('Amplitude: {0:0.4f}({1:0.4f})'.format(tup[0], perr_res[idx][0]))
    print('Center: {0:0.4f}({1:0.4f})'.format(tup[1], perr_res[idx][1]))
    print('Width: {0:0.4f}({1:0.4f})'.format(tup[2], perr_res[idx][2]))

print('\nOffset: {0:0.4f}({1:0.4f})'.format(popt[-1], perr[-1]))

plt.figure()
plt.plot(x, y, '-', linewidth=4, label='Data')

for idx, p in enumerate(popt[:-1].reshape(int((len(popt)-1)/3),3).tolist()):
    A, m, sig = p
    plt.plot(x, gaussian(x, A, m, sig) + popt[-1], linewidth=2, label='Peak{2}: ${0:0.3f} \, {1:0.1f}\,$Resolution'.format(m, m/(2.355*sig), idx))

plt.plot(x, multi_peaks(x, *popt), 'r--', linewidth=2, label='Fit')
plt.legend()
plt.show()

