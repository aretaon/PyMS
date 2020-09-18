# -*- coding: utf-8 -*-
"""
CalculateMW script to calculate peptide mass from amino acid sequence

@author: User
"""

_aa_masses = """A   71.03711   71.0788
R   156.10111   156.1875
N   114.04293   114.1038
D   115.02694   115.0886
C   103.00919   103.1388
E   129.04259   129.1155
Q   128.05858   128.1307
G   57.02146   57.0519
H   137.05891   137.1411
I   113.08406   113.1594
L   113.08406   113.1594
K   128.09496   128.1741
M   131.04049   131.1926
F   147.06841   147.1766
P   97.05276   97.1167
S   87.03203   87.0782
T   101.04768   101.1051
W   186.07931   186.2132
Y   163.06333   163.176
V   99.06841   99.1326"""

_label_names = ['15N',]

# masses from https://www.sciencedirect.com/science/article/abs/pii/0891391958900858
_label_masses =\
"""A   1.000036
R   4.000144
N   2.000072
D   1.000036
C   1.000036
E   1.000036
Q   2.000072
G   1.000036
H   3.000108
I   1.000036
L   1.000036
K   2.000072
M   1.000036
F   1.000036
P   1.000036
S   1.000036
T   1.000036
W   2.000072
Y   1.000036
V   1.000036"""

def CalculateMW(sequence, mono=False, returnAverage = False, labeling=None):
    """
    Calculate an accurate average or monoisotopic mass from an amino acid
    sequence

    args:
        sequence(str): The amino acid sequence
    kwargs:
        mono(bool): Return monoisotopic mass (default: False)
        returnAverage(bool): Return the mass divided by the sequence length (default: False)
        labeling(str): Return masses after metabolic labeling. Possible values are: None, '15N'. (Default: None)
    returns:
        mass(float): The calculated mass
    """

    ## Generate a modification dictionary if labeling is requested ##
    aminoacids = []
    mass = []
    if labeling is not None:
        try:
            label_idx = _label_names.index(labeling)
            aminoacids = []
            mass = []
            for line in _label_masses.split('\n'):
                line = line.split()
                aminoacids.append(line[0])
                mass.append(float(line[label_idx+1]))

            label_dict = dict(zip(aminoacids, mass))

        except:
            raise Exception('The label name provided was not found in the DB')

    ## Generate the main dictionary mapping amino acids names to masses ##
    aminoacids = []
    mass = []
    for line in _aa_masses.split('\n'):
        if mono == True:
            x, y, _ = line.split()
        else:
            x, _, y = line.split()
        aminoacids.append(x)
        if labeling is not None:
            y = float(y) + label_dict[x] # shift the mass by the label mass
        mass.append(float(y))
    use_dict = dict(zip(aminoacids, mass))

    # Sequentially sum up the masses of the amino acids
    mass = 18.01056 # condensation water
    for aa in sequence:
        if aa not in use_dict.keys():
            raise Exception('Amino acid {} not found in mass-list.'.format(aa))
        mass += use_dict[aa]

    if returnAverage is True:
        return mass/len(sequence)
    else:
        return mass