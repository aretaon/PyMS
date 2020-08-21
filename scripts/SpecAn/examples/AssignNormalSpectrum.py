# -*- coding: utf-8 -*-
"""
Annotate spectra with custom peptide masses
Created on Fri Oct 13 12:59:06 2017
@author: aretaon

"""
import os, sys
import matplotlib.pyplot as plt
sys.path.append(r'../..')
import SpecAn as sa

# peptide sequence to generate fragment ions from
peptide = 'GIGAVLKVLTTGLPALISWIKRKRQQ'
# precursor charge
prec_ch = 4
mods = [[-1,25]] # mods as mass, pos
# set the start for the mgf file
mgfPath = 'NormalSpectrum.mgf'
mgf_spec_num = 666
# error to tolerate during assignment of peaks
ppm = [-10, 10]


#%%
fig, ([ax1, ax2]) = plt.subplots(nrows=2,ncols=1,sharex=True)

# generate fragment ions from the peptide
ions2desc = sa.SpectrumAssignment.IonsFromSequence(peptide, # peptide
                               ['a', 'b', 'y'], # ion types
                               [1,4], # min, max charge
                               'monoisotopic', # mass type
                               2000, # max mass
                               mods=mods) # mods

# log all generated peptides to file
with open('peptides.log', 'w') as f:
    f.write('\t'.join(['Mass', 'Peptide', 'IonType', 'Length', 'Charge', 'PepType', 'Mods']) + '\n')
    for idx, i in enumerate(ions2desc.keys()):
        f.write('{}\t{}\n'.format(i, '\t'.join([str(i) for i in ions2desc[i]])))

# read the mgf file and extract an offset for each spectrum number (for faster retrieval)
with open(mgfPath, 'r') as f:
    spectrum2offset = sa.SpectrumReader.IndexMGF(f)
    
    # get the actual spectrum as list of lists mapping mz to intensity
    spectrum = sa.SpectrumReader.ReadSpectrum(mgf_spec_num, # spectrum
                                    f, # file handle
                                    spectrum2offset) # spectrum dict

Assignments = sa.SpectrumAssignment.AssignAndPlotPSM(spectrum, ions2desc, ppm, ax=ax1)

for idx, color in enumerate(Assignments['Color']):
    ax2.scatter(Assignments['AssignmentMz'][idx],
                Assignments['AssignmentError'][idx],
                color=color,
                alpha=0.5,
                s=10
                )
  
ax2.plot([100, 715],
         [0,0],
         linestyle = '--',
         lw=1,
         color='k'
         )

ax2.set_xlim([0,1.2*max(Assignments['AssignmentMz'])])
ax2.set_ylim(ppm)
ax2.set_ylabel('Error (ppm)')
ax2.set_xlabel('m/z')

plt.show()