# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 17:44:08 2017

@author: aretaon
"""

import pandas as pd
import os
import sys
import numpy as np

import matplotlib.pyplot as plt

def calc_rel_iBAQ(file, titles1, titles2, target_genes, normalise_to=None,
                  imputation = False):
    """
    Calculate the relative iBAQ value of a MaxQuant output.

    Keyword arguments:
    file -- the relative path to the input file
    titles1 -- a list of titles as given by MQ for the condition of interest
    titles2 -- same for the wildtype condition
    target_genes -- list of full or partly gene names to look for in the final
    analysis
    normalise_to -- gene name to normalise iBAQs to. NOrmalised to total iBAQ
    if None
    imputation -- impute zeros by normal distribution sampling (True or False)
    """

    # turn off warnings for chained dataframe assignments
    # this is raised by the append to dataframe operation in the last 2
    # lines of this function
    pd.options.mode.chained_assignment = None  # default='warn'

    # for pandas read_csv doc see
    # http://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html
    # returns a pandas dataframe object
    data = pd.read_csv(os.getcwd() + os.path.sep + file, sep='\t')

    # replace 0 by NaN
    data.replace(0, np.nan, inplace=True)
    # use the Gene names as index
    data.set_index(data['Gene names'], inplace=True)

    # avoid calculation a the inverse of an array of nans by checking first
    # if there are any reverse proteins in the data
    if not data['Reverse'].notnull().values.any():
        remain = data['Potential contaminant'] != '+'
    else:
        # remove contaminants and decoys
        # the truth statements must be grouped by parentheses, see
        # http://stackoverflow.com/questions/34531416/comparing-dtyped-float64-array-with-a-scalar-of-type-bool-in-pandas-datafram#34531543
        remain = (data['Reverse'] != '+') &\
            (data['Potential contaminant'] != '+')

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

    # print('===== Total iBAQs =====')
    # print(pex_proteins['iBAQ'])

    if imputation:

        # generate a copy of pex proteins to work on
        p = pex_proteins

        # initialise series to collect data
        mean_full = pd.Series()
        std_full = pd.Series()
        mean_subset = pd.Series()
        std_subset = pd.Series()

        # initialise the figure object
        fig, ax = plt.subplots(1, len(titles1 + titles2), sharex=True, sharey=True, squeeze=False)

        for idx, title in enumerate(titles1 + titles2):

            ax[0][idx].set_title(title)

            # take log2 from the filtered data containing the iBAQs of interest
            log2_iBAQs = np.log2(filtered["iBAQ %s" % title])

            log2_iBAQs.plot.hist(bins=100,
                                 range=(0, 40),
                                 ax=ax[0][idx],
                                 color='k',
                                 alpha=0.5)

            # calculate std and mean of original distribution
            std_full[title] = log2_iBAQs.std()
            mean_full[title] = log2_iBAQs.mean()

            # calculate the percentile for 1% of the data
            quan = log2_iBAQs.quantile(q=0.01)
            # define subset with less than 1% of the values
            subset_iBAQs = log2_iBAQs[log2_iBAQs < quan]

            subset_iBAQs.plot.hist(bins=100,
                                   range=(0,40),
                                   ax=ax[0][idx],
                                   color = 'r')

            # do some basic statistics to compute the distribution
            std_subset[title] = subset_iBAQs.std()
            mean_subset[title] = subset_iBAQs.mean()

            p["iBAQ %s" % title].fillna(2**np.random.normal(loc=mean_subset[title],
                                                            scale=std_subset[title]),
                                        inplace=True)

        d = {'mean_full': mean_full,
             'std_full': std_full,
             'mean_subset': mean_subset,
             'std_subset': std_subset}

        # create the dataframe
        ImputationInfo = pd.DataFrame(d)

        # save information
        ImputationInfo.to_csv('ImputationInfo.csv')
        plt.show()

    else:
        p = pex_proteins.fillna(0)


    # all relative iBAQ calculation is done in p while the results are returned
    # as part of pex_proteins.
    # Thus, the input data is returned extended by the relative_iBAQ column

    # manipulate the iBAQ values
    for index in range(len(titles1)):
        title1 = titles1[index]
        title2 = titles2[index]

        if normalise_to:

            # check if entry is available in data
            if not p[p['Gene names'] == normalise_to].empty:

                # dataframe query returns dataframe with one entry and NaNs
                # reduce to single value by summing up
                norm_iBAQ_t1 = p[p['Gene names'] == normalise_to]["iBAQ %s" % title1].sum()
                norm_iBAQ_t2 = p[p['Gene names'] == normalise_to]["iBAQ %s" % title2].sum()

                # calculate the relative iBAQ values
                rel_iBAQ_t1 = p["iBAQ %s" % title1] / norm_iBAQ_t1
                rel_iBAQ_t2 = p["iBAQ %s" % title2] / norm_iBAQ_t2

            else:
                sys.exit("You provided a gene name to normalise to but" +
                         "it doesnt exist in the data")

        else:
            # generate a total iBAQ by summing up the filtered iBAQs per title
            total_iBAQ_t1 = filtered['iBAQ %s' % title1].sum()
            total_iBAQ_t2 = filtered['iBAQ %s' % title2].sum()

            # calculate the relative iBAQ values
            rel_iBAQ_t1 = p["iBAQ %s" % title1] / total_iBAQ_t1
            rel_iBAQ_t2 = p["iBAQ %s" % title2] / total_iBAQ_t2

        # append to dataframe
        pex_proteins["rel_iBAQ %s" % title1] = rel_iBAQ_t1
        pex_proteins["rel_iBAQ %s" % title2] = rel_iBAQ_t2

    return pex_proteins

