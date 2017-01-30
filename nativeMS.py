# -*- coding: utf-8 -*-


def ExtractPeakrangeAsList(fname, start=1, end='end', returntype='tuple',
                           exclude=[(0, 0)], points='all'):
    """
    Extract a range of spectra from an mzML file and return them
    -----------

    Keyword arguments:

    fname:
        name of the input file \n

    start:
        the id of the first spectrum to return (counting starts at 1)

    end:
        the last spectrum to return, leave blank for last spectrum

    returntype:
        sets the export formatting.
        tuple returns a tuple with (all_mz, all_intensity)\n
        list returns a list of binary tuples [(mz, int), (mz, int), ...]\n
        tuple is required for digitizing

    exclude:
        takes a list of (start,end) tuples to be ignored during import.
        Useful for removing of spraying errors.\n

    points:
        takes an integer number of points to extract within the spectrum.
        Useful to reduce computation time

    Returns: mzlist, intlist
        mzlist: m/z coordinate\n
        intlist: intensity
    """

    import pymzml
    import numpy as np

    try:
        msrun = pymzml.run.Reader(fname)
    except:
        print('You must provide a valid filename')
        return

    if end == 'end':
        end = msrun.info['spectrum_count']
    else:
        end = int(end)

    if start == 'start':
        start = 1
    else:
        start = int(start)

    if exclude == [(0, 0)]:
        print('exporting from %(a)s to %(b)s' % {'a': start, 'b': end})
    else:
        print('exporting from %(a)s to %(b)s with exclusions' % {'a': start,
                                                                 'b': end})

    # Select a subset of all spectra to reduce computation time
    if points =='all':
        PointsInExtraction = int(end) - int(start)
    elif isinstance(points, int):
        if points < end:
            PointsInExtraction = points
        else:
            # avoid oversampling
            PointsInExtraction = int(end) - int(start)
    else:
        print('You must provide a valid integer for points\n'+
              'Using all data points')
        PointsInExtraction = int(end) - int(start)   

    SpectraToExtract = np.linspace(start,end,PointsInExtraction)
    SpectraToExtract = [int(round(r, 0)) for r in SpectraToExtract]

    if returntype == 'tuple':
        peakslist = []
        for i in SpectraToExtract:
            actual_peak = msrun[i].peaks
            peakslist.append(actual_peak)

        results = {}

        for peaks in peakslist:
            for key, value in peaks:
                results[key] = results.get(key, 0) + value

        results = list(results.items())

        return sorted(results)

    elif returntype == 'list':
        mzlist = []
        intlist = []
        for i in SpectraToExtract:
            for exclusion in exclude:
                if i not in range(exclusion[0], exclusion[1]):
                    for j in msrun[i].mz:
                        mzlist.append(j)
                    for k in msrun[i].i:
                        intlist.append(k)

        return mzlist, intlist
    else:
        print('valid export options are list or tuple!')
        return


def PlotPeaklist(peaklist, start='start', end='end', inputtype='tuple', plot='show'):
    """
    Returns a graphical representation of a spectrum list of tuples
    ------

    Keyword arguments:
        peaklist -- name of the variable containing the tuple\n
        inputtype -- either list or tuple
            tuple is a tuple with (all_mz, all_intensity)\n
            list is a list of binary tuples [(mz, int), (mz, int), ...]\n
        returned by ExtractPeakrangeAsList with the value export=list\n
        plot -- 'show' for direct plotting, return for return of the plot-object
        sg: True or False, whether filtered by savitzky golay filtering
    """

    import pylab as plt
    import sys

    if inputtype == 'list':
        intens = [b for a, b in peaklist]
        mz = [a for a, b in peaklist]

    elif inputtype == 'tuple':
        mz, intens = peaklist

    else:
        sys.exit('Inputtype is either list or tuple!')

    if start == 'start':
        start = 0

    if end == 'end':
        end = max(mz)

    intens_out = []
    mz_out = []

    # remove data outside set boundaries
    for i in range(len(mz)):
        if mz[i] > start:
            if mz[i] < end:
                mz_out.append(mz[i])
                intens_out.append(intens[i])

    f, ax = plt.subplots()
    ax.set_ylim([0, 1.1*max(intens_out)])
    ax.set_xlim([start, end])
    # see   http://stackoverflow.com/questions/10665163/
    #       draw-rectangle-add-patch-in-pylab-mode
    ax.plot(mz_out, intens_out, 'black')
    ax.set_xlabel('m/z [Th]')
    ax.set_ylabel('Rel. Intensity')

    if plot == 'show':
        plt.show()
    elif plot == 'return':
        return ax.get_figure()
    else:
        sys.exit('Plot is either show or return')

#    plt.plot(mz, intens)
#    plt.xlabel('m/z')
#    plt.ylabel('counts')
#    plt.savefig('%(a)s_extracted_from%(b)s_to_%(c)s_spectrum.eps' %
#                {'a': fname[:-5], 'b': start, 'c': end})
#    plt.clf()


def DigitizePeaklist(peaklist, resolution=0.5, debug=False,  sg = False,
                     inputtype='tuple',
                     returntype='tuple'):
    """
    Bin the list of peaks so that resolution decreases
    --------

    Keyword arguments:
        peaklist -- input peaklist as list of tuples \n
        resolution -- size of the bins, default 0.5 \n
        debug -- True activates verbose output
        sg -- apply savitzky golay filtering
        inputtype -- tuple or list
        returntype -- tuple is a outputs two lists of mz and intensities, list
        outputs one list of tuples of mz, int
    """
    import numpy as np
    import sys

    import savitzky_golay as sg


    if inputtype == 'tuple':
        mz_in, intens_in = zip(*peaklist)

        mz = np.array(mz_in)
        counts_double = np.array(intens_in)
    elif inputtype == 'list':
        mz = np.array(peaklist[0])
        counts_double = np.array(peaklist[1])
    else:
        sys.exit('You must provide a valid inputtype (list or tuple)')

    bins = np.arange(0, max(mz) + 1, resolution)

    inds = np.digitize(mz, bins)

    if debug:
        print(mz)
        print(bins)
        for n in range(round(mz.size/resolution)):
            print(bins[inds[n]-1], "<=", mz[n], "<", bins[inds[n]])

    mz_bins_double = []
    for m in range(mz.size):
        mz_bins_double.append(bins[inds[m]-1])

    mz_bins_double = np.array(mz_bins_double)

    mz_counts_tuple = {}
    for index, mass in enumerate(mz_bins_double):
        if mass in mz_counts_tuple.keys():
            mz_counts_tuple[mass] += counts_double[index]
        else:
            mz_counts_tuple[mass] = counts_double[index]

    mz_bins = []
    counts = []

    for element in sorted(list(mz_counts_tuple.items())):
        mz_bins.append(element[0])
        counts.append(element[1])

    if sg:
        counts = sg.savitzky_golay(counts, 49, 3)

    if returntype == 'tuple':
        return mz_bins, counts

    if returntype == 'list':
        return zip(mz_bins, counts)


def ExportPeaklistToUnidec(peaklist, outdir, inputtype='tuple'):

    """
    Writes a txt-file readable by UniDec for further spectra processing
    ------

    Keyword arguments:

        peaklist -- input peaklist, formatted as tuple
        (see ExtractPeakrangeAsList) \n
        outdir -- output directory
        istuple -- set to False allows to plot the list of binary tuples
    """

    import sys

    if inputtype == 'list':
        intens = [b for a, b in peaklist]
        mz = [a for a, b in peaklist]
    elif inputtype == 'tuple':
        mz, intens = zip(*peaklist)
    else:
        sys.exit('You must provide a valid input type')

    with open(outdir, 'w') as f:
        for index, mass in enumerate(mz):
            f.write('%(a).9f %(b).9f\n' % {'a': mass, 'b': intens[index]})
    f.close()


def DisplayTIC(fname, start=1, end='end', xaxis='time'):

    """
    Displays the TIC of a file and indicates two positions
    -----

    Keyword arguments:
        fname -- name of input file (mzML)\n
        start -- position of the first marker\n
        end -- position of the second marker\n
        xaxis -- format of the x-axis; either time (time) or spectrum # (id)
    """

    import pymzml
    import matplotlib.pyplot as plt

    msrun = pymzml.run.Reader(fname)

    if start == 'start':
        start = 1
    else:
        start = int(start)

    if xaxis == 'time':

        if end == 'end':
            # a random big number to fill the plot until the end
            # an explicit endtime cannot be calculated before iterating
            # over spectra
            end = 10000000
        else:
            end = int(end)

        time = []
        intensity = []
        for spectrum in msrun:
            try:
                time.append(spectrum['scan start time']*60)
                intensity.append(spectrum['total ion current'])
            except:
                continue

        print(start, end)

        f, ax = plt.subplots()
        ax.set_ylim([0, 1.1*max(intensity)])
        ax.set_xlim([0, max(time)])
        ax.add_patch(plt.Rectangle((start, 0), end-start, 1.1*max(intensity),
                                   fc='#336699',  # foreground colour
                                   ec='none'))  # edge colour
        ax.plot(time, intensity, 'black')
        ax.set_xlabel('time[sec]')
        ax.set_ylabel('intensity')

        f.tight_layout()
        return f

    if xaxis == 'id':

        if end == 'end':
            end = msrun.info['spectrum_count']
        else:
            end = int(end)

        scan_id = []
        intensity = []
        for spectrum in msrun:
            try:
                # removes entry 'TIC' on the last spectrum
                if isinstance(spectrum['id'], int):
                    scan_id.append(spectrum['id'])
                    intensity.append(spectrum['total ion current'])
            except:
                continue

        f, ax = plt.subplots()
        ax.set_ylim([0, 1.1*max(intensity)])
        ax.set_xlim([0, max(scan_id)])
        # see   http://stackoverflow.com/questions/10665163/
        #       draw-rectangle-add-patch-in-pylab-mode
        ax.add_patch(plt.Rectangle((start, 0), end-start, 1.1*max(intensity),
                                   fc='#336699',
                                   ec='none'))
        ax.plot(scan_id, intensity, 'black')
        ax.set_xlabel('scan id')
        ax.set_ylabel('intensity')

        f.tight_layout()
        return f

def DisplayXIC(fname,mass2follow):
    """
    under construction

    """

    import pymzml
    import pylab as plt

    try:
        msrun = pymzml.run.Reader(fname)
    except:
        print('You must provide a valid filename')
        return

    timeDependentIntensities = []
    for spectrum in msrun:
        if spectrum['ms level'] == 1:
            for time in spectrum:
                print(time)

