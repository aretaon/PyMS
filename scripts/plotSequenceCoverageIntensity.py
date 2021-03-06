import pandas as pd
import numpy as np
import os, sys

import matplotlib.pyplot as plt
import argparse

sys.path.append(r'C:\Users\User\Documents\03_software\python')
import PyMS.modules.Sequences as Seq

parser = argparse.ArgumentParser()
parser.add_argument('txt_path', help="Path to the MaxQuant txt folder")
parser.add_argument('target_name', help="Name of the target protein as stored in the fasta file")
parser.add_argument('fasta_path', help="Path to fasta file")
parser.add_argument('--use_index', dest='aaindex', action='store_true', help='Use amino acid indices for x axis')
parser.add_argument('--no_index', dest='aaindex', action='store_false', help='Use resnames for x axis')
parser.set_defaults(aaindex=False)
args = parser.parse_args()

infile = os.path.join(args.txt_path, 'evidence.txt')
targetProtein = args.target_name

fastafile = args.fasta_path

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
    
    ax.set_xlabel('Sequence position')
    if args.aaindex is True:
        plt.xticks(np.arange(0, len(thisSequence), 50))
    else:
        plt.xticks(np.arange(len(thisIntensity)),
               list(thisSequence),
               rotation=0)
    
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