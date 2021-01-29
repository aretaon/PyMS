#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Grubbs outlier test in Python
@author: aretaon
"""

import numpy as np
import scipy.stats as stats

def grubbs_test_masked(x, alpha, mask):
    """
    Perform 1D Grubbs' outlier test with a masked array'

    Parameters
    ----------
    x : np.array
        Data to evaluate.
    alpha : float
        alpha significance level for outlier detection.
    mask : numpy.array (bool)
        boolean array functioning as mask. Must be the same lenth as x.
        Values in x at the same position as True in mask will be ignored.

    Returns
    -------
    true_idx : int or None
        Integer of the detected outlier in the unmasked array. None if no outlier has been detected.

    """

    x = x[~mask]
    # print("Data as seen for Grubbs:\n{}".format(x))
    n = len(x)
    mean_x = np.mean(x)
    sd_x = np.std(x)
    numerator_idx = np.argmax(abs(x-mean_x))
    # print('Numerator idx = {}'.format(numerator_idx))
    numerator = max(abs(x-mean_x))
    g_calculated = numerator/sd_x
    print("Grubbs Calculated Value at N={}: {}".format(n, g_calculated))
    t_value = stats.t.ppf(1 - alpha / (2 * n), n - 2)
    g_critical = ((n - 1) * np.sqrt(np.square(t_value))) / (np.sqrt(n) * np.sqrt(n - 2 + np.square(t_value)))
    print("Grubbs Critical Value:",g_critical)
    if g_critical > g_calculated:
        print("Accept null hypothesis and conclude that there is no outliers")
        return None
    else:
        print("Reject null hypothesis and conclude that there is an outliers")

        masked_i = 0
        for true_idx, m in enumerate(mask):
            if masked_i == numerator_idx:
                break
            if m == False:
                masked_i += 1
        print("The outlier is {} at index position {}\n".format(x[numerator_idx], true_idx))

        return true_idx