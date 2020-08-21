def ReadFasta(fname):
    """
    Returns ID and sequence from a multi-item fasta-file
    
    Keyword arguments:
        fname -- path to fasta file
    
    Returns:
        accessions -- IDs of the proteins as list
        
        sequences -- Sequences of the proteins as list
    """
    with open(fname, 'r') as f:
        
        # initialise lists
        ID_list = []
        sequence_list = []
        
        # initialise sequence as boolean
        sequence = None
        
        # read all lines from file
        for line in f.readlines():
            line = line.strip()
            # the > is the identifier for a defline
            if line.startswith('>'):
                # use the inherent booleaness of string
                # every string is automatically a boolean TRUE
                if sequence:
                    sequence_list.append(sequence)  
                if '|' in line:
                    # take the part after the vertical line if present
                    ID = line.split('|')[1]
                else:
                    # else take the whole line without the tarting >
                    ID = line.strip()[1:]
                    
                # add the ID to the list of IDs
                ID_list.append(ID)    
                    
                # only initialise sequence if an ID has be found before
                sequence = ''
            else:
                # append the line to the sequence
                sequence += line
        
        sequence_list.append(sequence) # append the last sequence

        return ID_list, sequence_list



def peptide_cutter(sequence, enzyme='trypsin'):
    """
    Split a peptide string according to defined enzymatic cleavage rules.
    Regex for enzymes were taken from pyteomics
    (https://pyteomics.readthedocs.io/en/latest/api/parser.html#pyteomics.parser.expasy_rules)
    
    Keyword arguments:
        sequence(str) -- path to fasta file
        enzyme(str) -- enzyme to use, leave empty for full ist
    
    Returns:
        peptides(list) -- list of peptide sequence strings
    """
    import re
    
    enzymes = {  'arg-c': 'R',
                 'asp-n': '\\w(?=D)',
                 'bnps-skatole': 'W',
                 'caspase 1': '(?<=[FWYL]\\w[HAT])D(?=[^PEDQKR])',
                 'caspase 2': '(?<=DVA)D(?=[^PEDQKR])',
                 'caspase 3': '(?<=DMQ)D(?=[^PEDQKR])',
                 'caspase 4': '(?<=LEV)D(?=[^PEDQKR])',
                 'caspase 5': '(?<=[LW]EH)D',
                 'caspase 6': '(?<=VE[HI])D(?=[^PEDQKR])',
                 'caspase 7': '(?<=DEV)D(?=[^PEDQKR])',
                 'caspase 8': '(?<=[IL]ET)D(?=[^PEDQKR])',
                 'caspase 9': '(?<=LEH)D',
                 'caspase 10': '(?<=IEA)D',
                 'chymotrypsin high specificity': '([FY](?=[^P]))|(W(?=[^MP]))',
                 'chymotrypsin low specificity': '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
                 'clostripain': 'R',
                 'cnbr': 'M',
                 'enterokinase': '(?<=[DE]{3})K',
                 'factor xa': '(?<=[AFGILTVM][DE]G)R',
                 'formic acid': 'D',
                 'glutamyl endopeptidase': 'E',
                 'granzyme b': '(?<=IEP)D',
                 'hydroxylamine': 'N(?=G)',
                 'iodosobenzoic acid': 'W',
                 'lysc': 'K',
                 'ntcb': '\\w(?=C)',
                 'pepsin ph1.3': '((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|((?<=[^HKR][^P])[FLWY](?=\\w[^P]))',
                 'pepsin ph2.0': '((?<=[^HKR][^P])[^R](?=[FL][^P]))|((?<=[^HKR][^P])[FL](?=\\w[^P]))',
                 'proline endopeptidase': '(?<=[HKR])P(?=[^P])',
                 'proteinase k': '[AEFILTVWY]',
                 'staphylococcal peptidase i': '(?<=[^E])E',
                 'thermolysin': '[^DE](?=[AFILMV])',
                 'thrombin': '((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
                 'trypsin': '([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))'
                 }
    
    
    if enzyme not in enzymes.keys():
        raise Exception('Please use one of these enzymes: {}'.format('\n- '.join(enzymes.keys())))
        
    protein = re.sub(enzymes[enzyme],'\n', sequence)
    return protein.split()