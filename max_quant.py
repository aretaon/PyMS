# -*- coding: utf-8 -*-
"""
Methods related to MaxQuant analysis

@author: aretaon
"""

def calc_rel_iBAQ(file, titles1, titles2, target_genes):
    """
    Calculate the relative iBAQ value of a MaxQuant output.

    Keyword arguments:
    file -- the relative path to the input file
    titles1 -- a list of titles as given by MQ for the condition of interest
    titles2 -- same for the wildtype condition
    target_genes -- list of full or partly gene names to look for in the final
    analysis
    """

    # turn off warnings for chained dataframe assignments
    # this is raised by the append to dataframe operation in the last 2
    # lines of this function
    pd.options.mode.chained_assignment = None  # default='warn'

    # for pandas read_csv doc see
    # http://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html
    # returns a pandas dataframe object
    data = pd.read_csv(os.getcwd() + os.path.sep + file, sep='\t')

    # remove contaminants and decoys
    # the truth statements must be grouped by parentheses, see
    # http://stackoverflow.com/questions/34531416/comparing-dtyped-float64-array-with-a-scalar-of-type-bool-in-pandas-datafram#34531543
    remain = (data['Reverse'] != '+') & (data['Potential contaminant'] != '+')
    filtered = data[remain]

    try:
        len(titles1) == len(titles2)
    except:
        sys.exit("length of title lists must match!")

    # generate a gene names list to filter
    accepted = []
    for element in filtered['Gene names']:
        if isinstance(element, str):
            for gene_name in target_genes:
                if gene_name in element:
                    accepted.append(element)

    # reduce the filtered list to just the Pex-Proteins
    pex_proteins = filtered[filtered['Gene names'].isin(accepted)]

    for index in range(len(titles1)):
        title1 = titles1[index]
        title2 = titles2[index]

        # generate a total iBAQ by summing up the filtered iBAQs per title
        total_iBAQ_t1 = sum(filtered['iBAQ %s' % title1])
        total_iBAQ_t2 = sum(filtered['iBAQ %s' % title2])

        # calculate the relative iBAQ values
        rel_iBAQ_t1 = pex_proteins["iBAQ %s" % title1] / total_iBAQ_t1
        rel_iBAQ_t2 = pex_proteins["iBAQ %s" % title2] / total_iBAQ_t2

        # append to dataframe
        pex_proteins["rel_iBAQ %s" % title1] = rel_iBAQ_t1
        pex_proteins["rel_iBAQ %s" % title2] = rel_iBAQ_t2

    return pex_proteins
