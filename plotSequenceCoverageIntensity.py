import pandas as pd
import numpy as np
import sys

import matplotlib.pyplot as plt

sys.path.append(r'Z:\03_software\python')
from PyMS import Sequences as Seq

infile = 'evidence.txt'
targetProtein = 'His6-Syt1'

fastafile = 'His6_syt-1.fasta'

data = pd.read_csv(infile, delimiter='\t')
data = data[data['Proteins'] == targetProtein]
data = data[data['Intensity'].notnull()]

IDList, sequenceList = Seq.ReadFasta(fastafile)

thisSequence = sequenceList[IDList.index(targetProtein)]

intensities = np.zeros(len(thisSequence))

fig, ax = plt.subplots(1)
plt.xticks(np.arange(len(thisIntensity)),
           list(thisSequence),
           rotation=0)

for exp in set(data['Experiment']):
    expData = data[data['Experiment'] == exp]
        
    for row in expData.iterrows():
        peptide = row[1]['Sequence']
        start = thisSequence.index(peptide)
        end = start + len(peptide) - 1
        position = np.hstack((np.zeros(start),
                            np.ones(len(peptide)),
                            np.zeros(len(thisSequence) - end - 1)))
        thisIntensity = position * (row[1]['Intensity'])
        intensities += thisIntensity
    
    relIntensities = intensities / intensities.max()
    ax.plot(np.arange(len(thisIntensity)), relIntensities, label=exp)

ax.legend()
ax.set_xlabel('Sequence')
ax.set_ylabel('Rel Cumulative Intensity')
ax.set_ylim(bottom=0)
plt.show()