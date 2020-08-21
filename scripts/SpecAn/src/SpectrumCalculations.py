# -*- coding: utf-8 -*-
"""
Library of functions associated with ion mass calculation fram a sequence
Created on Wed Nov 22 09:26:33 2017

@author: User
"""

############################################
# Calculation of fragment peptide sequences
############################################

def NTermPeptides(parent, start=1, end=None):
    """
    return all n-terminal peptides from start (including)
    to end (non-including)
    """
    if end == None:
        end = len(parent)

    nterm = [parent[0:j] for j in range(start, end)]
    return nterm


def CTermPeptides(parent, start=1, end=None):
    """
    return all c-terminal peptides from start (non-including)
    to end (including)
    """
    if end == None:
        end = len(parent) 
    cterm = [parent[i:len(parent)] for i in range(start, end)]
    return cterm

def XPeptides(peptide1, xlink1, peptide2, xlink2, CTerm=False):
    """
    Calculate all fragmentation ion sequences for a cross-link peptide
    and return them sorted by whether they contain the cross-link or not
    
    :params: peptide1: sequence of the alpha peptide of the cross-link
    :params: xlink1: relative position of the cross-link within peptide1
    :params: peptide2: sequence of the beta peptide of the cross-link ion
    :params: xlink2:  relative position of the cross-link within peptide2
    :params: CTerm: Calc Cterminal peptides if true, else Nterminal
    """
    if CTerm:
        """
        Return all c-terminal peptides from a cross-linked sequence
        """
        # truncation of peptide 1
        unmod1 = CTermPeptides(peptide1, start=xlink1)
        mod1 = CTermPeptides(peptide1, end=xlink1)
        # return the indices of the beginning and the end of the substrings
        unmod1_pos = [xlink1+1, len(peptide1)]
        mod1_pos = [2, len(peptide1)]
        
        # truncation of peptide2
        unmod2 = CTermPeptides(peptide2, start=xlink2)
        mod2 = CTermPeptides(peptide2, end=xlink2)
        unmod2_pos = [xlink1+1, len(peptide2)]
        mod2_pos = [2, len(peptide2)]

    else:
        """
        Return all n-terminal peptides from a cross-linked sequence
        """
        # truncation of peptide1
        unmod1 = NTermPeptides(peptide1, end=xlink1)
        mod1 = NTermPeptides(peptide1, start=xlink1)
        # return the indices of the beginning and the end of the substrings
        unmod1_pos = [1, xlink1-1]
        mod1_pos = [1, len(peptide1)-1]
        
        # truncation of peptide2
        unmod2 = NTermPeptides(peptide2, end=xlink2)
        mod2 = NTermPeptides(peptide2, start=xlink2)
        unmod2_pos = [1, xlink2-1]
        mod2_pos = [xlink2, len(peptide2)]
        
    # return positions as single list of lists
    positions = [unmod1_pos,
                 mod1_pos,
                 unmod2_pos,
                 mod2_pos]
    
    return unmod1, mod1, unmod2, mod2, positions

############################################
# Calculation of masses of sequences
############################################

def BIonMass(peptide, aa_dict, charge, mods, peptide_type='alpha'):
    """
    Calculate the m/z of a b-ion with or without a modification
    
    Args:
        peptide: sequence of the peptide to generate mass from
        aa_dict: dict mapping amino acid residue names to masses
        charge:  charge of the peptide
        mods (optional): variable modificaitons of type [[mass1_pep1, pos1_pep1], [mass2_pep1, pos2_pep1]]
        peptide_type (optional):  either 'alpha' or 'beta'
    
    Returns:
        MassToCharge: mass of the y-ion
        Desc: description of the ion like [sequence, type, position,
              charge, peptide_type, Modifications
                     
    Raises:
        Exception: If peptide_type is specified other than 'alpha' or 'beta'
        
    """
    # check peptide_type
    if not peptide_type in ['alpha', 'beta']:
        raise Exception('You have to specify peptide_type as "alpha" or "beta"!')
    ion_string = 'B' if peptide_type == 'alpha' else 'b'
    
    # set initial values
    MassToCharge = 0
    moddesc = ''
    
    # calculate the total mass of the peptide    
    for idx, aa in enumerate(peptide):
        # add the modification mass if the iterator has passed the mod position
        for mod in mods:
            if idx + 1 == mod[1]:
                MassToCharge += mod[0]
                moddesc += str(mod[0]) + '({}) '.format(mod[1])
        MassToCharge += aa_dict[aa]
        
    # plus H20 from condensation - OH- loss from b-ion formation
    # proton charges = total charge - 1 charge from b-ion formation
    MassToCharge += 18.01528 - 17.00734 + (charge-1) * 1.00794
    MassToCharge /= charge # charging effect

    Desc = [peptide, ion_string, len(peptide), charge,
            peptide_type, moddesc]

    return MassToCharge, Desc

def AIonMass(peptide, aa_dict, charge, mods=[], peptide_type='alpha'):
    """
    Calculate the m/z of an a-ion with or without a modification
    
    Args:
        peptide: sequence of the peptide to generate mass from
        aa_dict: dict mapping amino acid residue names to masses
        charge:  charge of the peptide
        mods (optional): variable modificaitons of type [[mass1_pep1, pos1_pep1], [mass2_pep1, pos2_pep1]]
        peptide_type (optional):  either 'alpha' or 'beta'
    
    Returns:
        MassToCharge: mass of the y-ion
        Desc: description of the ion like [sequence, type, position,
              charge, peptide_type, Modifications]
                     
    Raises:
        Exception: If peptide_type is specified other than 'alpha' or 'beta'
        
    """
    # check peptide_type
    if not peptide_type in ['alpha', 'beta']:
        raise Exception('You have to specify peptide_type as "alpha" or "beta"!')
    ion_string = 'A' if peptide_type == 'alpha' else 'a'

    # set initial values
    MassToCharge = 0
    moddesc = ''

    # calculate the total mass of the peptide
    for idx, aa in enumerate(peptide):
        # add the modification mass if the iterator has passed the mod position
        for mod in mods:
            if idx + 1 == mod[1]:
                MassToCharge += mod[0]
                moddesc += str(mod[0]) + '({}) '.format(mod[1])
        MassToCharge += aa_dict[aa]
        
    MassToCharge += -28.0101 + charge * 1.00794 # M + nH+
    MassToCharge /= charge # charging effect

    Desc =  [peptide, ion_string, len(peptide), charge,
             peptide_type, moddesc]

    return MassToCharge, Desc

def YIonMass(peptide, aa_dict, charge, mods=[], peptide_len=None, peptide_type='alpha'):
    """
    Calculate the m/z of an y-ion with or without a modification
    
    Args:
        peptide: sequence of the peptide to generate mass from
        aa_dict: dict mapping amino acid residue names to masses
        charge:  charge of the peptide
        mods (optional): variable modificaitons of type [[mass1_pep1, pos1_pep1], [mass2_pep1, pos2_pep1]]
        peptide_len (optional): total length of peptide the modification position is calculcated from
        peptide_type (optional):  either 'alpha' or 'beta'
    
    Returns:
        MassToCharge: mass of the y-ion
        Desc: description of the ion like [sequence, type, position,
              charge, peptide_type, Modifications
                     
    Raises:
        Exception: If not and peptide_len are specified together
        Exception: If peptide_type is specified other than 'alpha' or 'beta'
        
    """

    # check peptide_type
    if not peptide_type in ['alpha', 'beta']:
        raise Exception('You have to specify peptide_type as "alpha" or "beta"!')
    ion_string = 'Y' if peptide_type == 'alpha' else 'y'

    # check modifications
    if (peptide_len == None and mods != []) or\
       (peptide_len != None and mods == []):
           raise Exception('If specifying mods you have to specify peptide length')

    # set initial values
    MassToCharge = 0
    moddesc = ''
    
    # calculate the total mass of the peptide
    for idx, aa in enumerate(peptide):
        # add the modification mass if the iterator has passed the mod position
        for mod in mods:
            if idx == peptide_len - mod[1]:
                MassToCharge += mod[0]
                moddesc += str(mod[0]) + '({}) '.format(mod[1])
        MassToCharge += aa_dict[aa]
        
    # plus H2O from condensation plus proton mass per charge
    MassToCharge += 18.01528 + charge * 1.00794 # M + nH+
    MassToCharge /= charge # charging effect

    Desc = [peptide, ion_string, len(peptide), charge,
            peptide_type, moddesc]

    return MassToCharge, Desc

def ParentionMass(peptide, aa_dict, charge, mods, peptide_type=None):
    ion_string = 'M+{}H'.format(charge)
    mass = 0
    moddesc = ''
    for idx, aa in enumerate(peptide):
        # add the modification mass if the iterator has passed the mod position
        for mod in mods:
            if idx + 1 == mod[1]:
                mass += mod[0]
                moddesc += str(mod[0]) + '({}) '.format(mod[1])
        mass += aa_dict[aa]
        
    mass += 18.01528 + charge * 1.00794 # M + nH+
    mass /= charge
    desc = [peptide, ion_string, '', charge,
            peptide_type, moddesc]
    return mass, desc

def XParentionMass(pep1_mass, pep2_mass, peptide1, xlink1, peptide2,
                    xlink2, charge, mods1, mods2, peptide_type=None):
    ion_string = 'M+{}H'.format(charge)
    mass = sum((pep1_mass, pep2_mass))
    moddesc = ''
    for mod in mods1:
        mass += mod[1]
        moddesc += str(mod[0]) + '({}) '.format(mod[1])
    for mod in mods2:
        mass += mod[1]
        moddesc += str(mod[0]) + '({}) '.format(mod[1])

    mass += charge * 1.00794 # M + nH+
    mass /= charge
    peptide = '-'.join([str(x) for x in [peptide1, xlink1, peptide2, xlink2]])
    desc = [peptide, ion_string, '', charge,
            peptide_type, moddesc]
    return mass, desc
