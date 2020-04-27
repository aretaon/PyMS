# -*- coding: utf-8 -*-
"""
Plot Massign or MassLynx derived spectrum lists


    Keyword arguments:

    SpectrumNames:
        file names of raw-spectra
    PlotStart:
        m/z value to start plotting
    PlotEnd:
        m/z value to end plot
    Normalise:
        normalise plots to max intensity
    Outname:
        name of the output file (default Massign_simulation.pdf)

Created on Wed Feb  1 10:37:22 2017

@author: aretaon
"""
import os
import pylab as plt

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--folder",
                    help="path to the the folder containing raw-spectra")
parser.add_argument("spectra",
                    help="file name of raw-spectrum",
                    nargs='+')
parser.add_argument("--plotstart",
                    type=int, help="m/z value to start plotting")
parser.add_argument("--plotend",
                    type=int, help="m/z value to end plotting")
parser.add_argument("--normalise",
                    type=bool, help="normalise plots to max intensity")
parser.add_argument("--outname",
                    help="name of the output file (default Massign_simulation.pdf")

args = parser.parse_args()

SpectraList = args.spectra

if args.folder:
    Folder = args.folder
else:
    Folder = os.getcwd()

# squeeze = False avoids compressing of arrays -> returns always 2D array
f, ax = plt.subplots(len(SpectraList), sharex=True, sharey=True, squeeze=False)

RawSpectra = []

for SpecName in SpectraList:

    with open(os.path.join(Folder, SpecName), 'r') as f:
        s = f.readlines()

        print("Reading: {}".format(SpecName))

        Specmz = []
        SpecInt = []

        for line in s:
            # split the data lines into array with mz and one with intensity
            Thismz = float(line.split('\t')[0].strip())
            ThisIntensity = float(line.split('\t')[1].strip())

            if (args.plotstart and args.plotstart <= Thismz) and \
               (args.plotend and args.plotend >= Thismz):

                    Specmz.append(Thismz)
                    SpecInt.append(ThisIntensity)
                    
        # normalise to max intensity
        if args.normalise:
            maxSpecInt = max(SpecInt)
            SpecInt = [100*i/maxSpecInt for i in SpecInt]

        ThisSpec = (Specmz[:], SpecInt[:], SpecName)
        # store mz and intensity as tuple in list
        RawSpectra.append(ThisSpec)


print("Found {} raw spectra:".format(len(RawSpectra)))
for i in RawSpectra:
    print(i[2])

for idx, Spectrum in enumerate(RawSpectra):

    ######## Get Start and End ###########

    if args.plotstart:
        PlotStart = args.plotstart
    else:
        PlotStart = 1

    print(PlotStart)

    if args.plotend:
        PlotEnd = args.plotend
    else:
        PlotEnd = max(Spectrum[0])

    print(PlotEnd)

    mz = Spectrum[0]
    Intensity = Spectrum[1]
    SpecName = Spectrum[2]

    print("Name: {}".format(Spectrum[2]))

    ######### Plot ##############

    print(ax[idx][0])

    # since squeeze == false we always get 2D arrays
    ax[idx][0].set_xlim([PlotStart, PlotEnd])

    ax[idx][0].plot(mz, Intensity, 'black', label=SpecName)

    if args.normalise:
        ax[idx][0].set_ylabel('% Intensity')
    else:
        ax[idx][0].set_ylabel('Intensity')
    ax[idx][0].legend(fontsize=8)

ax[idx][0].set_xlabel('m/z [Th]')
if args.outname:
    outpath = os.path.join(Folder, args.outname)
else:
    outpath = os.path.join(Folder, 'Massign_plot.pdf')

#plt.show()
plt.savefig(outpath)
