# -*- coding: utf-8 -*-
"""
Calculate mean-residue weighted CD spectra

@author: aretaon
"""

data_in = [
           [r'filename.txt', # infile
            r'baseline.txt', # baseline file
            600, #max_HT
            r'fastafile.fasta', # fastafile
            'P63041', #ID in fasta
            1.003, # concentration mg/mL
            0.01], # pathlength_cm
#           [r'second_filename.txt', # infile
#            None, # no baseline
#            600, #max_HT
#            r'fastafile.fasta', # fastafile
#            'synaptobrevin_49to96', #ID in fasta
#            1.0077, # concentration mg/mL
#            0.01], # pathlength_cm
]

plotrange= [195,250]

smoothing = [7,3]

### DONT CHANGE FROM HERE ###

import numpy as np
import matplotlib.pyplot as plt
import sys, os
from scipy.signal import savgol_filter
sys.path.append(r'C:\Softwaretools')
from PyMS.modules import Sequences as Seq
from PyMS.modules import CalculateMW as calcMW
from PyMS.modules.CDSpec import spectrumClass


# Perform MRW-calculation for all spectra
for idx, d in enumerate(data_in):
    
    infile, baseline, max_HT, fastafile, ID, concentration_mg_mL, pathlength_cm = d
    
    myspectrum = spectrumClass.spectrum(start_nm=plotrange[0], stop_nm=plotrange[1])
    myspectrum.load_jasco(infile)

    if baseline != None:
        mybaseline = spectrumClass.spectrum(start_nm=plotrange[0], stop_nm=plotrange[1])
        mybaseline.load_jasco(baseline)
        myspectrum.substract(mybaseline, inplace=True)
        
    myspectrum.filter_by_detector(maxV=max_HT)

    # Get the protein sequence from the fasta-file
    found = False
    IDList, sequenceList = Seq.ReadFasta(fastafile)
    for idx, identifier in enumerate(IDList):
        if identifier == ID:
            found = True
            protein_sequence = sequenceList[idx]
    if found == False:
        raise Exception('Could not find entry in DB')
    
    # Calculate the mean residue mass of the protein from the sequence    
    mean_residue_mass = calcMW.mass_from_sequence(protein_sequence, mono=False, returnAverage=False) / (len(protein_sequence) -1)
    
    
    myspectrum.get_MRWE(mean_residue_mass, pathlength_cm, concentration_mg_mL)
   
    with open(os.path.splitext(d[0])[0] + '_MRW.txt', 'w') as out:
        out.write('Wavelength[nm]\tMRW ellipticity[deg*cm^2*dmol^-1]\n')
        for idx, wl in enumerate(myspectrum.mrwe_deg_cm_2_dmol):
            out.write('{}\t{}\n'.format(wl, myspectrum.mrwe_deg_cm_2_dmol[idx]))
    #plt.plot(wavelength_nm, CD_mdeg, color='red')
    plt.plot(myspectrum.wavelength_nm, myspectrum.mrwe_deg_cm_2_dmol/10**4, label=os.path.basename(d[0]), linestyle='-')
# in-script smoothing
    plt.plot(myspectrum.wavelength_nm[~np.isnan(myspectrum.mrwe_deg_cm_2_dmol)],
                           savgol_filter(myspectrum.mrwe_deg_cm_2_dmol[~np.isnan(myspectrum.mrwe_deg_cm_2_dmol)]/10000,
                                                smoothing[0],
                                                smoothing[1]),
                           label=os.path.basename(d[0]) + '_{},{}_savgol_smooth'.format(smoothing[0], smoothing[1]),
                           linestyle='-')

# plot a zero-line
plt.plot(np.arange(*plotrange), np.zeros(len(np.arange(*plotrange))), color='gray', linestyle='--', linewidth=1)

# polishing
plt.ylabel('$\Theta_{MRW} * 10^{-4}$ [deg $cm^{2}$ $dmol^{-1}$]')
plt.xlabel('Wavelength [nm]')
plt.legend()
plt.show()