#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Rayleigh limit charge calculator for native mass spectrometry

Returns the charge of a spherical droplet at Rayleigh limit based on the
molecular mass of a protein (assuming average protein density and spherical
protein dimensions)

Created on Sun Jan 29 20:29:55 2017

@author: areaton
"""

import numpy as np
import scipy.constants as spc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("mass", type=float, help="Mass (g/mol) for the Rayleigh charge to be returned")

args = parser.parse_args() # in g/mol

molar_mass = args.mass

surf_tension = 0.073 # in N / m (for water)
permittivity = 8.854187817*10**-12 # in C² / N m² (for vacuum)

density = 1410 + 145 * np.exp(-molar_mass/13) # calculated after
# Average protein density is a molecular-weight-dependent function
# HANNES FISCHER, IGOR POLIKARPOV AND ALDO F. CRAIEVICH

molecules_per_droplet = 1 # assuming ion evaporation model

molecule_mass = molar_mass / spc.Avogadro # mass for one molecule in g

molecule_mass = molecule_mass / 1000 # mass for one molecule in kg

volume = molecule_mass / density # in m³

radius = (3/(4*np.pi) * volume) ** (1/3) # in m

print("Radius of droplet is {:.2f} nm".format(radius*10**9))

rayleigh_limit = 8 * np.pi * (permittivity * surf_tension * radius ** 3) ** 0.5

print("Charge if droplet at Rayleigh limit: {} C".format(rayleigh_limit))

# 1 Coulomb is equivalent to the charge of approximately 
# 6.242×10^18 (1.036×10^−5 mol) protons
no_protons = 6.242 * 10 ** 18 * rayleigh_limit

print("No of protons at Rayleigh limit: {:.2f}".format(no_protons))