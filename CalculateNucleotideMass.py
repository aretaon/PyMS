# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:57:40 2019

@author: User
"""

# MW from http://biotools.nubic.northwestern.edu/OligoCalc.html#helpOD

nucleotide2mass = {'RNA': {'A': 329.21, 'U': 306.17, 'C': 305.18, 'G': 345.21},
                   'DNA': {'A': 313.21, 'T': 304.2, 'C': 289.18, 'G': 329.21}}

fiveprime_mod2mass = {'OH': -62,
                      'Monophosphate': 0,
                      'Triphosphate': 159}

threeprime_mod2mass = {'None': 0}

with open('nucleotide_modifications.txt', 'r') as inf:
    for line in inf:
        if not line.startswith('#'):
            line = [x.strip() for x in line.split('\t')]
            # is fiveprime
            if line[2] == 'X':
                fiveprime_mod2mass[line[0]] = float(line[1])
            # is threeprime
            if line[3] == 'X':
                threeprime_mod2mass[line[0]] = float(line[1])
##%
#
#sequence= 'GAUGAGGGGUCACAUAACCACACUGUUGACUACCUCAUUCCUGGUUUUUGGCCUCCACAUC'
#seq_type = 'RNA'
#fiveprime_mod = 'Cyanine (CY) 5'
#threeprime_mod = "3'-BHQ-1 (Black Hole Quencher)"

sequence= 'UGAGGUAGUAGGUUGUGUGGUU'
seq_type = 'RNA'
fiveprime_mod = 'OH'
threeprime_mod = 'None'

mass = 0
mass += fiveprime_mod2mass[fiveprime_mod]

for n in sequence:
    try:
        mass += nucleotide2mass[seq_type][n]
    except:
        raise Exception('Could not find nucleotide {}'.format(n))

mass += threeprime_mod2mass[threeprime_mod]

print("The MW of the {}-sequence {} with the 5'mod {} and 3'mod {} is: {}".format(seq_type, sequence, fiveprime_mod, threeprime_mod, mass))