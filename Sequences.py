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

