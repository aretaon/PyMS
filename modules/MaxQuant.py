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

def calc_rel_iBAQ(file, titles, target_gene_names, normalise_to=None,
                  imputation = False):
    """
    Calculate the relative iBAQ value of a MaxQuant output.

    Keyword arguments:
    file -- the relative path to the input file
    titles -- a list of expepriment names as given to MQ
    target_gene_names -- list of full or partly Protein IDs to look for in the final
    analysis, set to None for all genes
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
    # use the Protein IDs as index
    data.set_index(data['Protein IDs'], inplace=True)

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

    # Only return a subset of genes of interest if variable is set
    # else calculate relativie iBAQs for the whole dataset
    if target_gene_names != None:
        # generate a Protein IDs list to filter
        accepted = []
        for element in filtered['Protein IDs']:
            if isinstance(element, str):
                for gene_name in target_gene_names:
                    if gene_name in element:
                        accepted.append(element)
    
        # reduce the filtered list to just the Pex-Proteins
        trgt_genes = filtered[filtered['Protein IDs'].isin(accepted)]
    else:
        trgt_genes = filtered

    if imputation:

        # generate a copy of pex proteins to work on
        t = trgt_genes

        # initialise series to collect data
        mean_full = pd.Series()
        std_full = pd.Series()
        mean_subset = pd.Series()
        std_subset = pd.Series()

        # initialise the figure object
        fig, ax = plt.subplots(1, len(titles), sharex=True, sharey=True, squeeze=False)

        for idx, title in enumerate(titles):

            ax[0][idx].set_title(title)
            ax[0][idx].set_xlabel('Log2 iBAQ')
            

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

            # replace NaN by small numbers sampled from the quantile
            t["iBAQ %s" % title].fillna(2**np.random.normal(loc=mean_subset[title],
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
        plt.savefig('ImputationInfo.pdf')

    else:
        t = trgt_genes.fillna(0)


    # all relative iBAQ calculation is done in p while the results are returned
    # as part of trgt_genes.
    # Thus, the input data is returned extended by the relative_iBAQ column

    # manipulate the iBAQ values
    for title in titles:

        if normalise_to:

            # check if entry is available in data
            if not t[t['Protein IDs'] == normalise_to].empty:

                # dataframe query returns dataframe with one entry and NaNs
                # reduce to single value by summing up
                norm_iBAQ_t1 = t[t['Protein IDs'] == normalise_to]["iBAQ %s" % title].sum()

                # calculate the relative iBAQ values
                rel_iBAQ_t1 = t["iBAQ %s" % title] / norm_iBAQ_t1

            else:
                sys.exit("You provided a gene name to normalise to but" +
                         "it doesnt exist in the data")

        else:
            # generate a total iBAQ by summing up the filtered iBAQs per title
            total_iBAQ_t1 = filtered['iBAQ %s' % title].sum()

            # calculate the relative iBAQ values
            rel_iBAQ_t1 = t["iBAQ %s" % title] / total_iBAQ_t1

        # append to dataframe
        trgt_genes["rel_iBAQ_%s" % title] = rel_iBAQ_t1

    return trgt_genes

def ReadProteinGroups(file):
    """
    Open a MaxQuant protein groups file and return it as pandas dataframe
    """
    pg = pd.read_csv(file,
                   delimiter='\t')

    # avoid calculation a the inverse of an array of nans by checking first
    # if there are any reverse proteins in the data
    if not pg['Reverse'].notnull().values.any():
        remain = pg['Potential contaminant'] != '+'
    else:
        # remove contaminants and decoys
        # the truth statements must be grouped by parentheses, see
        # http://stackoverflow.com/questions/34531416/comparing-dtyped-float64-array-with-a-scalar-of-type-bool-in-pandas-datafram#34531543
        remain = (pg['Reverse'] != '+') &\
            (pg['Potential contaminant'] != '+')
    
    data = pg[remain].copy()
    
    # replace 0 by NaN
    data.replace(0, np.nan, inplace=True)

    return data

    