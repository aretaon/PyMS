#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 17:06:13 2017

@author: aretaon
"""

import itertools
import pandas as pd

stuff = []

subunits = ['SubA',
            'SubB',
            'SubC',
            'SubD',
            'SubE']
masses = [2, 3, 4, 5, 6]

subunits2masses = dict(zip(subunits, masses))

max_subunits = 10
max_mass = 30

pot_complexes = list(itertools.product(subunits, repeat=max_subunits))

def test_mass(comp):
    mass = 0
    for sub in comp:
        mass += subunits2masses[sub]
        if mass >= max_mass:
            return False
    return mass

def nice(comp):
    comp_string = ''
    for sub in subunits:
        comp_string += sub
        comp_string += '_{} '.format(comp.count(sub))
    return comp_string

true_complexes = []
true_masses = []

for comp in pot_complexes:
    mass = test_mass(comp)
    if mass:
        true_complexes.append(comp)
        true_masses.append(mass)

nice_complexes = []

for element in true_complexes:
        nice_complexes.append(nice(element))

data = pd.DataFrame.from_dict({'Complex' : nice_complexes,
                               'Mass' : true_masses})

data.to_excel('Subunit combinations.xlsx', index=False)