import pandas as pd
import numpy as np
import sys, os

import matplotlib.pyplot as plt

sys.path.append(r'Z:\03_software\python')
from PyMS import Sequences as Seq

os.chdir(r'C:\Users\User\Documents\02_experiments\07_instrument_tests\jb07b_native_MS_membrane_proteins\20191021_mq_analysis')

infile = 'evidence.txt'
targetProtein = 'AcrB_calixar'

fastafile = '20191021_AcrB.fasta'

data = pd.read_csv(infile, delimiter='\t')
data = data[data['Proteins'].str.contains(targetProtein)]
data = data[data['Intensity'].notnull()]

IDList, sequenceList = Seq.ReadFasta(fastafile)

thisSequence = sequenceList[IDList.index(targetProtein)]

def generate_plot(df, thisSequence, label=None):
    
    fig, ax = plt.subplots(1)
    intensities = np.zeros(len(thisSequence))

    for row in df.iterrows():
        peptide = row[1]['Sequence']
        start = thisSequence.index(peptide)
        end = start + len(peptide) - 1
        position = np.hstack((np.zeros(start),
                            np.ones(len(peptide)),
                            np.zeros(len(thisSequence) - end - 1)))
        thisIntensity = position * (row[1]['Intensity'])
        intensities += thisIntensity
    
    relIntensities = intensities / intensities.max()
    colors = ['red' if x <= 0 else 'green' for x in relIntensities]
    if label != None:
        ax.plot(np.arange(len(thisIntensity)), relIntensities, label=label)
    else:
        ax.plot(np.arange(len(thisIntensity)), relIntensities)
    for pos in np.arange(len(thisIntensity)):
        ax.plot(pos, relIntensities[pos], 'o', color=colors[pos], alpha=0.5)
        
    ax.legend()
#    ax.set_xlabel('Sequence')
#    plt.xticks(np.arange(len(thisIntensity)),
#            list(thisSequence),
#            rotation=0)
    
    ax.set_xlabel('Sequence position')
    plt.xticks(np.arange(0, len(thisSequence), 50))
    
    ax.set_ylabel('Rel Cumulative Intensity')
    ax.set_ylim(bottom=0)

# If experiments are defined, the plot will be generated for every experiment separately
if 'Experiment' in data.columns:
    for exp in set(data['Experiment']):
        expData = data[data['Experiment'] == exp]
        generate_plot(expData, thisSequence, exp)
# else all peptide information for the single portein will be plotted together
else:
    generate_plot(data, thisSequence)

plt.show()