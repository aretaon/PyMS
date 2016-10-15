# -*- coding: utf-8 -*-
"""
Methods related to statistics

@author: aretaon
"""

def basic_stats(data):
    """
    Return mean, standard deviation and standard error of the mean for pandas
    data series

    Keyword arguments:
    data: a pandas data series

    Returns:
    mean, std, SEM, N
    """
    print(data)
    mean = float(data.mean())
    std = float(data.std())
    # count the filled entries by subtracting the empty one from the size
    N = len(data) - data.isnull().sum()
    print(N)
    SEM = float(std/np.sqrt(N))

    return mean, std, SEM, N
