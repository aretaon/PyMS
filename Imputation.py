def lowImputer(a, percentile):
    """
    Imputes missing values (NaN) by sampling from the
    lowest percentile of values from the series
    """
    
    import numpy as np
    
    a = np.array(a)
    
    if sum(np.isnan(a)) == len(a):
        return a
    
    # calculate std and mean of original distribution
    std_full = np.nanstd(a)
    mean_full = np.nanmean(a)
    # calculate the percentile for 1% of the data
    perc = np.nanpercentile(a, percentile)

    print('a: %s' % a)
    
    b = a[~np.isnan(a)]
    
    subset = b[b < perc]
    # integer is required to reduce computational effort for distribution calc
    subset = subset.astype(int)
    
    print('subset: %s' % subset)
    
    # do some basic statistics to compute the distribution
    std_subset = subset.std()
    mean_subset = subset.mean()

    print('std_subs: %s' % std_subset)
    print('mean_subs: %s' % mean_subset)
    
    #Find indicies that you need to replace
    inds = np.where(np.isnan(a))

    # replace NaN by small numbers sampled from the quantile
    # a.fillna(2**np.random.normal(loc=mean_subset, scale=std_subset),
    #                             inplace=True)
            
    if std_subset == 0:
        a[inds] = mean_subset
    else:
        #Place column means in the indices. Align the arrays using take
        a[inds] = np.random.normal(loc=mean_subset, scale=std_subset)
    
    return a