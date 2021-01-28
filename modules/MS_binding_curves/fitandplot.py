# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 20:16:26 2020

@author: aretaon
"""

import os, glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from outliertests import grubbs_test_masked

#%%%% Gather data %%%%

os.chdir(r'/home/aretaon/mlu_cloud/Til_MSc/High resolution plot/Synapt')

# min and max x value of the plot
plt_range = [0,200]

# read the correction terms from csv and generate a dictionary mapping lipid names to
# correction factors
correction = pd.read_csv(r'./concentration_correction.csv', encoding="latin-1", usecols=[0,1], index_col=0).iloc[:,0].to_dict()
    
# collect data and assemble into dataframe object
all_data = [] 
for file in glob.glob('*/rep*/extracts.txt'): # use extracts.txt (contains extracted intensities)
    print(file)
    data = pd.read_csv(file, delimiter=' ', header=None)
    try:
        # if there is a x_axis.txt in the folder, read it
        x = np.genfromtxt(os.path.join(os.path.split(file)[0], 'x_axis.txt'))
    except:
        # if not, just continue (code usually breaks in the next line)
        continue
    # assign the lipid name from the filename
    data['lipid'], data['rep'] = os.path.split(file)[0].split(os.path.sep)
    # set the x values
    data['x'] = x * correction[os.path.split(file)[0].split(os.path.sep)[0]]
    # attach to the data stack
    all_data.append(data)
    
# concatenate all stacked yx-sets
df = pd.concat(all_data, ignore_index=False, axis=0, sort=False)
# reset the index to sequential integers
df = df.reset_index()

# make the single datasets accessible from the stack by grouping according
# to lipid name
grouped = df.groupby('lipid')

# add some defined colors for the lipids
color_dict = {'DOPS': 'orange',
              'DOPC': 'rosybrown',
              'DOPG': 'crimson',
              'DOPI': 'skyblue',
              'PI(4)P': 'royalblue',
              'PIP2': 'blueviolet',
              'PIP3': 'black'}

#%%%% Iterate through lipid groups %%%%

# this is the fit function
def hyperbola(x, kd):
#    return 1-x/(x+kd)
    # a*x**b
    return x/(x+kd)

for datacols, plotname in [((1,2,3), 'One Lipid Bound'), ((2,3), 'Two Lipids Bound'), ([3], 'Three Lipids Bound')]:
    
    # init the figure object
    fig, ax = plt.subplots(1,1,
                         figsize=(10, 10)) # determines the size of the figure in inches
    
    # Single binding model fitted to all bound states
    # sequentially plot all the data groups onto the figure object
    for idx, key in enumerate(grouped.groups.keys()):
    
        #%%%% Remove Outliers with Grubbs %%%%    
    
        # selecta subset of data belonging to one lipid
        d = grouped.get_group(key).fillna(0) # NaN are replaced with zero to allow summation
        # plot data points
        # use the second (i.e. number 1) column for plotting
        # change the number here for
        # 0: lipidfree
        # 1: first lipid
        # 2: 2nd lipid
        
        # sum all single peaks as in every peak there is at least one molecule binding
        # y_data = d[1] + d[2] + d[3]
        y_data = sum([d[i] for i in datacols])
    
        print('\n\n===============\Fitting {}\n===============\n'.format(key))
        
        # if key == 'PIP3':
        #     masked = d.x > 25
        #     masked = masked.to_numpy() # needs numps for proper indexing
        # else:
        masked = np.zeros(len(d.x), dtype=bool)
    
        popts = []
    
        while True:
            # fit the function to the datapoints and collect the optimised parameters
            # i.e. the result of the fit and the covariance of the optimised function
            # with regard to the datapoints
            popt, pcov = curve_fit(hyperbola, d.x[~masked], y_data[~masked], maxfev = 10000)
            popts.append(*popt)
        
            # calculate the deviation of all the points to the fit
            distance_to_fit = np.abs(y_data - hyperbola(d.x, *popt)).to_numpy()
            # print("Absolute distance to fit:\n{}".format(distance_to_fit))
            # print("Current mask:\n{}".format(masked))
            
            # obtain an outliner based on the masked array
            outlier_idx = grubbs_test_masked(distance_to_fit, 0.05, masked)
            
            if outlier_idx != None:
                masked[outlier_idx] = True # index here is array index beginning at 0
            else:
                break    
    
        #%%%% Fit and plot %%%%
    
        # refit the function with the reduced dataset
        popt, pcov = curve_fit(hyperbola, d.x[~masked], y_data[~masked], maxfev = 10000)
        popts.append(*popt)
    
        ax.scatter(d.x[~masked],
                    y_data[~masked],
                    color=color_dict[key])
        
        ax.plot(d.x[masked],
                y_data[masked],
                color=color_dict[key],
                fillstyle='none',
                linestyle='none',
                marker='o')
    
        kd = popt[0]
        err = np.sqrt(np.diag(pcov))
    
        # plot the optimised function as line onto the canvas
        ax.plot(np.linspace(*plt_range), hyperbola(np.linspace(*plt_range), kd), color=color_dict[key],
                 label='{0: <5}'.format(key) + r' $K_{app}=$' + ' {:10.2f}'.format(*popt) + r' $\mu M$')
        
        if len(popts)>1:
            lower_border = np.amax(np.array([hyperbola(np.linspace(*plt_range), kd) for kd in popts]), axis=0)
            upper_border = np.amin(np.array([hyperbola(np.linspace(*plt_range), kd) for kd in popts]), axis=0)
            
            ax.fill_between(np.linspace(*plt_range),
                             lower_border,
                             upper_border,
                             color=color_dict[key],
                             alpha=0.2)
            
        # add a legend with the lipid names to the top
        # plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=3,
        #         borderaxespad=0, frameon=False)
        # add a legend with the lipid names to the right
        plt.legend(loc='upper left', bbox_to_anchor= (1.1, 1.0), ncol=1,
                borderaxespad=0, frameon=False)
    
    #%%%% plotting options %%%%
    
    ax.set_title(plotname)
    # these are the axis labels
    ax.set_xlabel('Lipid concentration [$\mu$M]')
    ax.set_ylabel('Relative Intensity')
    ax.set_xlim(plt_range)
    ax.set_ylim([0,1])
    # ax.semilogx()

# makes plots nicer
plt.tight_layout()
# shows the plot
plt.show()

# plt.savefig('lipid_plot.png', dpi=400)
# plt.savefig('lipid_plot.pdf', dpi=400)
