# -*- coding: utf-8 -*-
"""
Plot Massign Simulations
-----------

    Keyword arguments:

    SimFolder:
        path to the the folder containing simulations and raw-spectrum
    SimNames:
        list of file names for simulations to visualise
    Spectrum:
        file name of raw-spectrum
    PlotStart:
        m/z value to start plotting
    PlotEnd:
        m/z value to end plot
    EnvelopeResolution:
        for how many points should the envelope be calculated
    Outname:
        name of the output file (default Massign_simulation.pdf)


Created on Tue Jan 31 15:43:26 2017

@author: aretaon
"""

import os
import pylab as plt
import re
import numpy as np

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--simfolder",
                    help="path to the the folder containing simulations and raw-spectrum")
parser.add_argument("--simnames",
                    help="list of file names for simulations to visualise",
                    nargs='+')
parser.add_argument("spectrum",
                    help="file name of raw-spectrum")
parser.add_argument("--plotstart",
                    type=int, help="m/z value to start plotting")
parser.add_argument("--plotend",
                    type=int, help="m/z value to end plotting")
parser.add_argument("--normalise",
                    type=bool, help="normalise plots to max intensity")        
parser.add_argument("--outname",
                    help="name of the output file (default Massign_simulation.pdf")


args = parser.parse_args()

if args.simfolder:
    SimFolder = args.simfolder
else:
    SimFolder = os.getcwd()

if args.simnames:
    SimNames = args.simnames

Spectrum = args.spectrum

EnvelopeResolution = 1500

f, ax = plt.subplots(len(SimNames), sharex=True, sharey=True)

######### Get data for raw-spectrum ##########

with open(os.path.join(SimFolder, Spectrum), 'r') as f:
    Spectrum = f.readlines()

Specmz = []
SpecInt = []

for line in Spectrum:
    # split the data lines into array with mz and one with intensity
    Thismz = float(line.split('\t')[0])
    ThisIntensity = float(line.split('\t')[1])

    if (args.plotstart and args.plotstart <= Thismz) and \
       (args.plotend and args.plotend >= Thismz):

            Specmz.append(Thismz)
            SpecInt.append(ThisIntensity)
        
# normalise to max intensity
if args.normalise:
    SpecInt = [100*i/max(SpecInt) for i in SpecInt]

########### Get data for every simulation and plot ###########

if args.plotstart:
    PlotStart = args.plotstart
else:
    PlotStart = 1

print("Plotstart: {}".format(PlotStart))

if args.plotend:
    PlotEnd = args.plotend
else:
    PlotEnd = max(Specmz)

print("Plotend: {}".format(PlotEnd))

for idx, SimName in enumerate(SimNames):

    SimPath = os.path.join(SimFolder, SimName)

    with open(SimPath, 'r') as f:
        Simulation = f.readlines()

    mz = []
    Intensity = []
    Header = ''

    FindInHeader = re.compile(r'mass: (\d*) \+\/\- (\d*).*chargestates : (\d*) to (\d*) max (\d*).*FWHM : (\d*).*attachments : (\d*).*amplitude : ([-+]?[0-9]*\.?[0-9]+).*center : ([-+]?[0-9]*\.?[0-9]+).*standard deviation : ([-+]?[0-9]*\.?[0-9]+)')

    for line in Simulation:
        # sort the first lines not starting with a number
        if not line[0].isnumeric():
            Header = Header + line.strip()
        else:
            # split the data lines into array with mz and one with intensity
            Thismz = float(line.split('\t')[0].strip())
            ThisIntensity = float(line.split('\t')[1].strip())
            
            if (args.plotstart and args.plotstart <= Thismz) and \
               (args.plotend and args.plotend >= Thismz):

                    mz.append(Thismz)
                    Intensity.append(ThisIntensity)
                    
    # normalise to max intensity
    if args.normalise:
        Intensity = [100*i/max(Intensity) for i in Intensity]

    ######### Extract Information from Header ###########

    HeaderInfo = FindInHeader.match(Header)

    # only assign variables if header is found
    if HeaderInfo:

        Mass = float(HeaderInfo.group(1))
        MassError = float(HeaderInfo.group(2))
        ChargeMin = HeaderInfo.group(3)
        ChargeMax = HeaderInfo.group(4)
        FWHM = HeaderInfo.group(6)
        Attachments = HeaderInfo.group(7)
        EnvAmplitude = float(HeaderInfo.group(8))
        EnvCenter = float(HeaderInfo.group(9))
        EnvStd = float(HeaderInfo.group(10))

    ########## Generate envelope ############

    def gaussian(x, amp, mu, sig):
        return amp*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

    print('Amplitude: {}'.format(EnvAmplitude))
    print('Std: {}'.format(EnvStd))
    print('Center: {}'.format(EnvCenter))

    Envelope = gaussian(np.linspace(PlotStart, PlotEnd, EnvelopeResolution),
                        EnvAmplitude,
                        EnvCenter,
                        200*EnvStd)

    ######### Plot ##############

    ax[idx].set_xlim([PlotStart, PlotEnd])

    ax[idx].plot(Specmz, SpecInt, 'black')
    
    label = "Mass: {} +/- {} kDa".format(round(Mass/1000,2),
                                         round(MassError/1000, 2))
    ax[idx].plot(mz, Intensity, 'red', label=label)
    #ax[idx].plot(np.linspace(PlotStart, PlotEnd, EnvelopeResolution), Envelope)

    if args.normalise:
        ax[idx].set_ylabel('% Intensity')
    else:
        ax[idx].set_ylabel('Intensity')
        
    # only label simulation if header is found
    if HeaderInfo:
        ax[idx].legend(fontsize=8)

ax[idx].set_xlabel('m/z [Th]')
if args.outname:
    outpath = os.path.join(SimFolder, args.outname)
else:
    outpath = os.path.join(SimFolder, 'Massign_simulation.pdf')

plt.savefig(outpath)