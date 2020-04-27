# -*- coding: utf-8 -*-
"""
CalculateMW script to calculate peptide mass from amino acid sequence

@author: User
"""

aa_masses = """A	71.03711	71.0788
R	156.10111	156.1875
N	114.04293	114.1038
D	115.02694	115.0886
C	103.00919	103.1388
E	129.04259	129.1155
Q	128.05858	128.1307
G	57.02146	57.0519
H	137.05891	137.1411
I	113.08406	113.1594
L	113.08406	113.1594
K	128.09496	128.1741
M	131.04049	131.1926
F	147.06841	147.1766
P	97.05276	97.1167
S	87.03203	87.0782
T	101.04768	101.1051
W	186.07931	186.2132
Y	163.06333	163.176
V	99.06841	99.1326"""

aminoacids = []
mono_mass = []
average_mass = []

for line in aa_masses.split('\n'):
    x, y, z = line.split('\t')
    aminoacids.append(x)
    mono_mass.append(float(y))
    average_mass.append(float(z))

aa2monomass = dict(zip(aminoacids, mono_mass))
aa2averagemass = dict(zip(aminoacids, average_mass))

def mass_from_sequence(sequence, mono=False, returnAverage = False):
    sequence = sequence.strip()
    if mono == True:
        use_dict = aa2monomass
    else:
        use_dict = aa2averagemass
    mass = 18 # condensation water
    for aa in sequence:
        if aa not in use_dict.keys():
            raise Exception('Amino acid {} not found in mass-list.'.format(aa))
        mass += use_dict[aa]
        
    if returnAverage is True:
        return mass/len(sequence)
    else:    
        return mass