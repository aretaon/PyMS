# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 09:33:18 2016

@author: aretaon
"""

import argparse

import os

import nativeMS as nMS
import manipulate_text as mt

import matplotlib.pyplot as plt

import time

import csv

# to export single figures into pdf
from matplotlib.backends.backend_pdf import PdfPages

thisdate = time.strftime("%Y%m%d")

#################### Options ############################

ColsToPlot = ['Kommentar',
              'Puffer',
              'Probenname',
              'Collision Energy',
              'Transfer Collision Energy',
              'Trap Collision Energy',
              'Sample prep']

#################### Argparser ##########################

parser = argparse.ArgumentParser()

parser.add_argument("rawfiles", help="path to tsv-file to print spectra from",
                                type=str)

parser.add_argument('-m', '--mzmldir', help='Path dor dir storing mzml versions\
                                             of the files of interest. \n For\
                                             conversion unse the proteowizard\
                                             package.')

parser.add_argument('-t', '--title', help='Title for output pdf-file')

# generate parser object -> extract arguments from input
args = parser.parse_args()

# Read arguments from input
ParameterDataPath = args.rawfiles.strip()
print('Reading rawfiles from {}'.format(ParameterDataPath))

if args.mzmldir == None:
    mzmlPath = '/home/aretaon/AGW/massspec/_people/jbender/1_data/mzml-files'
else:
    mzmlPath = args.mzmldir.strip()

if args.title != None:
    main_title = args.title
else:
    main_title = '{}_extracted_spectra'.format(thisdate)

############## Extract Rawfile Names ####################

# open the excel file and extract the Rawfiles column
RawfilesNames = []
AnnText = []
with open(ParameterDataPath, 'r') as r:
    ParameterData = csv.DictReader(r,
                                   dialect='excel-tab')
    for row in ParameterData:
        RawfilesNames.append(row['Rawfile'])

        # extract the annotation information
        AnnTextTemp = []
        for Col in ColsToPlot:
            if Col in row.keys():
                AnnTextTemp.append('{}: {}'.format(Col, row[Col]))

        AnnText.append(AnnTextTemp)
# extract rawfile names from the first column of the excel file
RawfilesNames = [r.split('.')[0] for r in RawfilesNames]

print('Extracted rawfiles: \n{}'.format('\n'.join(RawfilesNames)))

############# Check if corresponding mzml exist ##########
# List of mzml files in the mzmlPath to extract spectra from
mzmlNames = []
# iterate through rawnames and check for existence of mzml files
for rawfile in RawfilesNames:
    found = False
    # check different version of the mzml extension
    PossibleRawnames = [rawfile + ext for ext in ['.mzml', '.mzML']]
    for r in PossibleRawnames:
        if found is True:
            continue  # avoid using both extensions
        FullPath = os.path.join(mzmlPath, r)
        if os.path.isfile(FullPath):
            mzmlNames.append(FullPath)
            found = True
    if found is not True:
        print('Didnt find a mzMl file for {} at {}'.format(rawfile, FullPath))
        mzmlNames.append(None)

############################ Plotting ####################

# List of figure objects for plotting
PlotList = []

# fancier plots
plt.style.use('ggplot')

# save the pdf in the current working dir using main_title as title
pp = PdfPages(os.path.join(os.getcwd(), mt.alphanum_string(main_title)) + '.pdf')

for idx, m in enumerate(mzmlNames):

    # print only spectra for rawfiles w/ mzml file
    if m is not None:

        print('Parsing {}'.format(os.path.basename(m)))

        # extract the raw Peaklist from the mzml file
        PeakList = nMS.ExtractPeakrangeAsList(fname=m,
                                              start=1,
                                              end='end',
                                              returntype='tuple',
                                              points=100,
                                              exclude=[(0, 0)])

        # digitize and smoothe (savitzky golay) peaklist to clear the view
        PeakList = nMS.DigitizePeaklist(PeakList,
                                        resolution=0.5,
                                        debug=False,
                                        sg=True,
                                        inputtype='tuple',
                                        returntype='tuple')

        # fetch a figure object for the spectrum plot
        fig = nMS.PlotPeaklist(PeakList,
                               start='start',
                               end='end',
                               inputtype='tuple',
                               plot='return')

        # get current axis from figure
        ax = fig.gca()

        # give the figure a title
        ax.set_title(os.path.basename(m))

        # add annotation to the plots in the top right corner of the plot
        plt.text(1, 1, '\n'.join(AnnText[idx]),
                 va='top',
                 ha='right',
                 fontsize=6,
                 transform=ax.transAxes,
                 bbox=dict(facecolor='red',
                           alpha=0.5))

        # save plots to pdf
        pp.savefig(fig)

    # catch rawfile w/o any mzml file == if m is None
    else:
        
        # Set empty variable to print no spectrum
        plt.clf() # clear figure
        
        # give the figure a title
        ax.set_title(os.path.basename(RawfilesNames[idx]))

        # write error message instead of plot to pdf file
        plt.text(0.5, 0.5, 'mzML file not found',
                 va='top',
                 ha='right',
                 fontsize=12,
                 transform=ax.transAxes)

        # add annotation to the plots in the top right corner of the plot
        plt.text(1, 1, '\n'.join(AnnText[idx]),
                 va='top',
                 ha='right',
                 fontsize=6,
                 transform=ax.transAxes,
                 bbox=dict(facecolor='red',
                           alpha=0.5))
        
        # save plots to pdf
        pp.savefig(fig)

# close handler for pdf
pp.close()