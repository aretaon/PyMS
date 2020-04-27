def lookupPDBChains(pdbPath):
    """
    Return all chain IDs present in a PDB file
    
    Args:
        pdbPath: path to PDB file
    Returns:
        chains: list of chain identifiers
    """

    allChains = set()

    with open(pdbPath) as f:
        for line in f.readlines():
            if line.startswith('COMPND'):
                if 'CHAIN:' in line:
                    IDString = line[line.find('CHAIN:') + 6:]
                    IDString = IDString.strip().strip(';').split(',')
                    IDs = set([x.strip() for x in IDString])
                    allChains = allChains | set(IDs)
    
    print(sorted(list(allChains)))
    
    return sorted(list(allChains))