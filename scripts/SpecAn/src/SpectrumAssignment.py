# -*- coding: utf-8 -*-
"""
Annotate spectra with results from cross-link identification searches.
Part of the CroCo project.
Created on Fri Oct 13 12:59:06 2017
@author: aretaon

"""
import os

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import csv
import numpy as np

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from . import SpectrumCalculations as Calc
from . import SpectrumReader as Reader
from . import XlinkPic

# set the start for relative paths if script run in standalone
binPath = '../../bin'


#%% Start of IonsFromSequence
def IonsFromSequence(sequence, ion_types, charge, mass_type, max_mass, mods=[]):
    """
    Calculates theoretical ions for an amino acid
    sequence with modifications.

    Args:
        sequence: Peptide sequence
        ion_types: list of which ions to return (a/b/y)
        charge state of the returned ions
        mass_type: mono or average
        max_mass: maximal mass to calculate ions from (e.g. max quad mass)
        mods (optional): modifications to consider (of type [[mass1_pep1, pos1_pep1], [mass2_pep1, pos2_pep1]])
    
    Returns:
        ions2desc: Dictionary mapping ion masses to their corresponding description in the form [sequence, type, position,
              charge, peptide_type, Modifications]
    
    Raises:
        FileNotFoundError: If the amino acid dict is not found
    """

    # table contains amino acid masses - H2O for better calculation
    aaPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                          './config/aa_masses.txt')
    try:
        with open(aaPath, mode='r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            if mass_type == 'monoisotopic':
                aa_dict = {rows[0]:float(rows[1]) for rows in reader}
            elif mass_type == 'average':
                aa_dict = {rows[0]:float(rows[2]) for rows in reader}
            else:
                print('No mass_type specified, assuming monoisotopic masses')
                aa_dict = {rows[0]:float(rows[1]) for rows in reader}
    except:
        raise FileNotFoundError('aa_masses.txt not found')

    all_ions = []
    all_desc = []

    for pep in Calc.NTermPeptides(sequence):
        if 'a' in ion_types:
            for c in range(charge[0], charge[-1]+1):
                mass, desc = Calc.AIonMass(pep, aa_dict, c, mods)
                all_ions.append(mass)
                all_desc.append(desc)
        if 'b' in ion_types:
            for c in range(charge[0], charge[-1]+1):
                mass, desc = Calc.BIonMass(pep, aa_dict, c, mods)
                all_ions.append(mass)
                all_desc.append(desc)

    for pep in Calc.CTermPeptides(sequence):
        if 'y' in ion_types:
            for c in range(charge[0], charge[-1]+1):
                mass, desc = Calc.YIonMass(pep, aa_dict, c, mods=mods, peptide_len=(len(sequence)))
                # print(sequence)
                # print(len(sequence))
                # print(desc)
                all_ions.append(mass)
                all_desc.append(desc)

    for c in range(charge[0], charge[1]+1):
        mass, desc = Calc.ParentionMass(sequence, aa_dict, c, mods)
        all_ions.append(mass)
        all_desc.append(desc)

    ions2desc = dict(zip(all_ions, all_desc))

    return ions2desc

#%% Start of IonsFromXlinkSequence
def IonsFromXlinkSequence(peptide1, xlink1, peptide2, xlink2, xlinker_mod, ion_types, charge, mass_type, max_mass, mods1=[], mods2=[]):
    """
    Calculates theoretical ions for an amino acid
    sequence with modifications.

    Args:
        peptide1: Peptide1 sequence
        xlink1: Relative cross-link position in peptide1
        peptide2: Peptide2 sequence
        xlink2: Relative cross-link possition peptide 2
        xlinker_mod: Modification mass of the crosslinker
        ion_types: list of which ions to return (a/b/y)
        charge state of the returned ions
        mass_type: monoisotopic or average
        mods1: variable modifications to consider in peptide1
        mods2: variable modifications to consider in pepetide2
        max_mass: maximum allowed mass (e.g. max m/z of quad)
        
    Returns:
        ions2desc: Dict mapping ion mz to description in the form [sequence, type, position,
              charge, peptide_type, Modifications]
    """

    #%% CalcABIons
    def CalcABions(peptides, peptide_type, ion_types, charge, aa_dict,
                   xpepdesc, xpepmass, mods):
        """
        Calculate theoretical masses and corresponding descriptions
        of a-ions for set of peptides sequences. Append the masses and descriptions
        to all_ions and all_dessc if mass < max_mass

        Args:
            peptides: List of peptide sequences
            peptide_type: string "alpha" or "beta" (first or last peptide truncated respectively)
            ion_types: list of ion types to consider (e.g. ['b', 'y'])
            charge: list of [min_charge, max_charge]
            aa_dict: dict linking one-character ID to mass
            xpepdesc: description to add to the xlinked peptides
            xpepmass: mass of the cross-linker plus second peptide
            mods: [[[mass1_1, mass2_1], [pos1_1, pos2_1]]
        """
        for pep in peptides:
            if 'a' in ion_types:
                for c in range(charge[0], charge[-1]+1):
                    mass, desc = Calc.AIonMass(pep, aa_dict, c,
                                                mods, peptide_type)
                    if mass < max_mass:
                        all_ions.append(mass + xpepmass/c) # treat second peptide as constant adduct
                        desc[0] += xpepdesc
                        all_desc.append(desc)

            if 'b' in ion_types:
                for c in range(charge[0], charge[-1]+1):
                    mass, desc = Calc.BIonMass(pep, aa_dict, c,
                                                mods, peptide_type)
                    if mass < max_mass:
                        all_ions.append(mass + xpepmass/c) # treat second peptide as constant adduct
                        desc[0] += xpepdesc
                        all_desc.append(desc)


    #%% CalcYions
    def CalcYions(peptides, peptide_type, ion_types, charge, aa_dict,
                   xpepdesc, xpepmass, mods=[], peptide_len=None):
        """
        Calculate theoretical masses and corresponding descriptions
        of y-ions for set of peptides sequences. Append the masses and descriptions
        to all_ions and all_dessc if mass < max_mass

        Args:
            peptides: List of peptide sequences
            peptide_type: string "alpha" or "beta" (first or last peptide truncated respectively)
            ion_types: list of ion types to consider (e.g. ['b', 'y'])
            charge: list of [min_charge, max_charge]
            aa_dict: dict linking one-character ID to mass
            xpepdesc: description to add to the xlinked peptides
            xpepmass: mass of the cross-linker plus second peptide
            mods: [[[mass1_1, mass2_1], [pos1_1, pos2_1]]
        """

        for pep in peptides:
            if 'y' in ion_types:
                for c in range(charge[0], charge[-1]+1):
                    mass, desc = Calc.YIonMass(pep, aa_dict, c,
                                               mods=mods,
                                               peptide_len=peptide_len,
                                               peptide_type=peptide_type)
                    if mass < max_mass:
                        all_ions.append(mass + xpepmass/c)
                        desc[0] += xpepdesc
                        all_desc.append(desc)

    #%% Continue IonsFromXlinkSequence

    # table contains amino acid masses - H2O for better calculation
    aaPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                          './config/aa_masses.txt')
    with open(aaPath, mode='r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        if mass_type == 'monoisotopic':
            aa_dict = {rows[0]:float(rows[1]) for rows in reader}
        elif mass_type == 'average':
            aa_dict = {rows[0]:float(rows[2]) for rows in reader}
        else:
            print('No or wrong mass_type specified, assuming monoisotopic masses')

    all_ions = []
    all_desc = []

    # collect sequences from modified and unmodified peptides (n-term of CID)
    nterm_unmod1, nterm_mod1, nterm_unmod2, nterm_mod2, nterm_positions=\
        Calc.XPeptides(peptide1, xlink1, peptide2, xlink2, CTerm=False)
    # C-terminal of CID break
    cterm_unmod1, cterm_mod1, cterm_unmod2, cterm_mod2, cterm_positions =\
        Calc.XPeptides(peptide1, xlink1, peptide2, xlink2, CTerm=True)

    print('Unmodified nterm peptides: {}'.format((', '.join(nterm_unmod1 +
                                                            nterm_unmod2))))

    print('Unmodified cterm peptides: {}'.format((', '.join(cterm_unmod1 +
                                                            cterm_unmod2))))

    print('Modified nterm peptides: {}'.format((', '.join(nterm_mod1 +
                                                            nterm_mod2))))

    print('Modified cterm peptides: {}'.format((', '.join(cterm_mod1 +
                                                            cterm_mod2))))

    ######################################################
    # Nterminal non xlinked peptides
    ######################################################

    nterm_unmod1 = Calc.NTermPeptides(peptide1, end=xlink1)
    nterm_unmod2 = Calc.NTermPeptides(peptide2, end=xlink2)

    CalcABions(nterm_unmod1, 'alpha', ion_types, charge, aa_dict,
               '', # xpepdesc
               0, # xpepmass
               mods1)
    CalcABions(nterm_unmod2, 'beta', ion_types, charge, aa_dict,
               '', # xpepdesc
               0, # xpepmass
               mods2)

    ######################################################
    # Nterminal xlinked peptides
    ######################################################

    # calculate mass of modification with peptide2 and xlinker
    pep2_mass = sum([aa_dict[aa] for aa in peptide2]) + 18.01528
    xpepmass1 = xlinker_mod + pep2_mass
    print('Mass of modification for peptide 1: {}'.format(xpepmass1))
    # same with peptide1
    pep1_mass = sum([aa_dict[aa] for aa in peptide1]) + 18.01528
    xpepmass2 = xlinker_mod + pep1_mass
    print('Mass of modification for peptide 2: {}'.format(xpepmass2))
    # define descriptions for peptide2 bound to (truncated) peptide1
    xpepdesc1 = '-{}-{}-{}'.format(xlink1, peptide2, xlink2)
    # same for peptide1 bound to (truncated) peptide2
    xpepdesc2 = '-{}-{}-{}'.format(xlink2, peptide1, xlink1)

    CalcABions(nterm_mod1, 'alpha', ion_types, charge, aa_dict,
               xpepdesc1, # xpepdesc
               xpepmass1, # xpepmass
               mods1)

    CalcABions(nterm_mod2, 'beta', ion_types, charge, aa_dict,
               xpepdesc2, # xpepdesc
               xpepmass2, # xpepmass
               mods2)

    ######################################################
    # Cterminal non-xlinked peptides
    ######################################################

    CalcYions(cterm_unmod1, 'alpha', ion_types, charge, aa_dict,
               '', # xpepdesc
               0, # xpepmass
               mods=mods1,
               peptide_len=len(peptide1))

    CalcYions(cterm_unmod2, 'beta', ion_types, charge, aa_dict,
               '', # xpepdesc
               0, # xpepmass
               mods=mods2,
               peptide_len=len(peptide2))

    ######################################################
    # Cterminal xlinked peptides
    ######################################################

    CalcYions(cterm_mod1, 'alpha', ion_types, charge, aa_dict,
               xpepdesc1, # xpepdesc
               xpepmass1, # xpepmass
               mods=mods1,
               peptide_len=len(peptide1))

    CalcYions(cterm_mod2, 'beta', ion_types, charge, aa_dict,
               xpepdesc2, # xpepdesc
               xpepmass2, # xpepmass
               mods=mods2,
               peptide_len=len(peptide2))

    ######################################################
    # Parent ions
    ######################################################

    for c in range(charge[0], charge[1]+1):
        mass, desc = Calc.XParentionMass(pep1_mass, pep2_mass,
                                          peptide1, xlink1,\
                                          peptide2, xlink2, c, mods1, mods2)
        if mass < max_mass:
            all_ions.append(mass)
            all_desc.append(desc)

    ions2desc = dict(zip(all_ions, all_desc))

    return ions2desc


#%%

def AssignAndPlotPSM(spectrum, ions2desc, ppm, ax=None, max_mz=None):
    """
    Annotate a given spectrum with a given sequence of theoretical
    mz and descriptions.

    Args:
        spectrum: list of lists of experimental mz and intensity
        ion2desc: dict mapping theoretical mz to description
        ppm: error allowed during assign in parts-per-million
        axis (optional): axis to plot the figure on, default currrent axis

    Returns:
        Assignments: Dict with descriptions of assignments as keys and lists
                     of the corresponding values in order as values
    """

    mz_list = spectrum[0]
    intensity_list = spectrum[1]
    ion_list = sorted(ions2desc.keys())

    if ax is None:
        ax = plt.gca()

    last = 0

    def ColorFromDesc(desc_list):
        """
        Return the plink color scheme for a given ion description
        """
        if desc_list[1] == 'A':
            return 'yellow'
        elif desc_list[1] == 'a':
            return 'brown'
        elif desc_list[1] == 'B':
            return 'green'
        elif desc_list[1] == 'b':
            return 'orange'
        elif desc_list[1] == 'Y':
            return 'red'
        elif desc_list[1] == 'y':
            return 'purple'
        elif 'M' in  desc_list[1]:
            return 'cyan'

    AssignmentError = []
    AssignmentMz = []
    Color = []
    Descriptions = []
    
    for idx, mz in enumerate(mz_list):
        # store axis-object in variable for return
        ax.plot([mz, mz], # x1, x2
                [0, intensity_list[idx]], # y1, y2
                'k-', # style
                lw=0.5 # line width
                )
        for idx, theor_mz in enumerate(ion_list[last:]):
            if theor_mz >= mz * (1+ppm[0]/10**6):
                if theor_mz <= mz * (1+ppm[1]/10**6):
                    # add the relative error of the assignments to the list
                    AssignmentError.append((mz-theor_mz)/mz * 10**6)
                    AssignmentMz.append(theor_mz)
                    Descriptions.append(ions2desc[theor_mz])
                    Color.append(ColorFromDesc(ions2desc[theor_mz]))

                    ax.plot([theor_mz, theor_mz],
                               [0, intensity_list[idx]],
                               ColorFromDesc(ions2desc[theor_mz]),
                               lw=1,
                               alpha=0.5)
                                        
                    ax.text(theor_mz,
                              intensity_list[idx]+1000,
                              '${0}_{{{1}}}^{{{2}+}} {3}$'.format(ions2desc[theor_mz][1],
                                                              ions2desc[theor_mz][2],
                                                              ions2desc[theor_mz][3],
                                                              ions2desc[theor_mz][5]),
                              {'ha': 'left', 'va': 'bottom'},
                              rotation=90,
                              color=ColorFromDesc(ions2desc[theor_mz]))
                    last += idx
    if max_mz == None:
        ax.set_xlim(min(mz_list),max(mz_list))
    else:
        try:
            ax.set_xlim(min(mz_list),max_mz)
        except:
            print('Couldnt set xlim to {}.\nUsing highest mz measured insted'.format(max_mz))
            ax.set_xlim(min(mz_list),max(mz_list))
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')

    if Descriptions != []:
        [Sequences, Types, Positions, Charges, PepTypes, Modifications] =\
            list(zip(*Descriptions))

        Assignments = {'AssignmentError': AssignmentError,
                       'AssignmentMz': AssignmentMz,
                       'Color': Color,
                       'Sequences': Sequences,
                       'Types': Types,
                       'Positions': Positions,
                       'Charges': Charges,
                       'PepTypes': PepTypes,
                       'Modifications': Modifications}
        
        return Assignments
    else:
        raise Exception('Could not assign any peak with the given specifications')


def PlotHist(data, limits, ax=None):
    """
    Plots a histogram with a fitted gaussion distribution over a set
    of data on a given axis

    Args:
        data: a list of values to plot a histogram from
        limits: range of the plot as list [lower, upper]
        ax: axis to plot the data on
    """
    if ax is None:
        ax = plt.gca()
    n, bins, patches = ax.hist(assignment_error,
                         50, # bins
                         normed=1,
                         alpha=0.8)
    mean = np.mean(assignment_error)
    variance = np.var(assignment_error)
    sigma = np.sqrt(variance)

    # add a 'best fit' line
    y = mlab.normpdf(bins, mean, sigma)
    ax.plot(bins, y, 'r--', linewidth=1)
    ax.axvline(mean, color='r', linewidth=1)

    ax.set_xlabel('Error in ppm')
    ax.set_ylabel('Probability Density')
    ax.set_xlim(limits)

#%%

if __name__ == '__main__':

    fig, ([ax1, ax2], [ax3, ax4]) = plt.subplots(nrows=2,
                                                 ncols=2,
                                                 sharex='col',
                                                 gridspec_kw = {'height_ratios':[3, 1],
                                                                'width_ratios': [3,1]})

    ax4.axis('off')

    ppm = [-20, 20]

    assignment_error = []

    peptide1 = 'KAWGNNQDGVVASQPAR'
    peptide2 =  'LKSSDAYK'
    xlink1 = 1
    xlink2 = 2
    prec_ch = 3

    ions2desc =  IonsFromXlinkSequence(peptide1,
                                       xlink1,
                                       peptide2,
                                       xlink2,
                                       138.068, #mass of xlinker
                                       ['a', 'b', 'y'], # ion types
                                       [1,3], # min, max charge
                                       'monoisotopic', # mass type
                                       2000, # max mass
                                       #[[5, 3], ], # mods1
                                       '', #mods1
                                       '') #mods2

    # ions2desc =  IonsFromSequence('PEPTIDE', # peptide
    #                               ['a', 'b', 'y'], # ion types
    #                               [1,1], # min, max charge
    #                               'monoisotopic', # mass type
    #                               2000, # max mass
    #                               [[300, 6]]) # mods

    with open('peptides.log', 'w') as f:
        f.write('\t'.join(['Mass', 'Peptide', 'IonType', 'Length', 'Charge', 'PepType', 'Mods']) + '\n')
        for idx, i in enumerate(ions2desc.keys()):
            f.write('{}\t{}\n'.format(i, '\t'.join([str(i) for i in ions2desc[i]])))

    mgfPath = os.path.join(binPath,
                           '../testdata/SV_plink/2017_08_04_SVs_BS3_16.mgf')
    with open(mgfPath, 'r') as f:
        spectrum2offset = Reader.IndexMGF(f)

        mz2intens = Reader.ReadSpectrum(18452, # spectrum
                                        f, # file handle
                                        spectrum2offset) # spectrum dict

    Assignments = AssignAndPlotPSM(mz2intens, ions2desc, ppm, ax=ax1)

    XlinkPic.PlotXlinkMap(peptide1,
                          peptide2,
                          xlink1,
                          xlink2,
                          prec_ch,
                          Assignments['Positions'],
                          Assignments['Types'],
                          Assignments['Charges'],
                          Assignments['PepTypes'],
                          fontsize=8,
                          ax=ax2)

    for idx, color in enumerate(Assignments['Color']):
        ax3.scatter(Assignments['AssignmentMz'][idx],
                    Assignments['AssignmentError'][idx],
                    color=color,
                    alpha=0.5,
                    s=10
                    )
  
    ax3.plot([0, 2000],
             [0,0],
             linestyle = '--',
             lw=1,
             color='k'
             )

    ax3.set_xlim([0,1.2*max(Assignments['AssignmentMz'])])
    ax3.set_ylim(ppm)
    ax3.set_ylabel('Error (ppm)')
    ax3.set_xlabel('m/z')

    #PlotHist(assignment_error, ppm, ax=ax[1])

    plt.show()