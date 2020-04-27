# -*- coding: utf-8 -*-
"""
Waters Metadata Plotter

Given a folder path the programme collects metadata from every Waters rawfile
living in this directory and generates an xls file from it

Created on Tue Mar  5 13:38:49 2019

@author: aretaon
"""


import argparse
import os
import re
import codecs
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", default=os.getcwd(), help="Path to the directory containing the raw files")
parser.add_argument("-o", "--outpath", default=os.getcwd(), help="Path to a folder to write the output")
parser.add_argument('-a', "--all", action="store_true", default=False, help="Print all columns and not only those that differ between the samples")
args = parser.parse_args()

## FOR TESTING UNCOMMENT HERE
#
#args.dir = r'C:\Users\User\Documents\03_software\python\PyMS\watersMetadataPlotter_testdata'
#args.all = False
#args.outpath = r'C:\Users\User\Documents\03_software\python\PyMS\watersMetadataPlotter_testdata'

##########

## Find the rawfiles to look at
rawfiles = list()

for root, dirs, files in os.walk(args.dir):
    for f in dirs:
        if f.endswith('.raw'):
            rawfiles.append(os.path.join(root, f))

base_rawfiles = [os.path.split(f)[1][:-4] for f in rawfiles]

## set the outfile
outname = '{:s}_complete_tune.csv'.format(os.path.basename(args.dir))

## parse _extern.inf and _HEADER.TXT for key value pairs
headerPattern = re.compile(r'\$\$ ([^:]+)\:(.*)')
data = list()
all_colnames = list()

def add2dict(d, k, v):
    if k in d.keys():
        d[k+'_2'] = v
    else:
        d[k] = v

for r in rawfiles:

    ParamDict = {}

    ### Collection data ###
    ParamDict['Rawfile'] = os.path.split(r)[1]
    all_colnames.append('Rawfile')

    # encoding has to be set to latin1 as MassLynx does not use UTF-8
    try:
        inhandler = codecs.open(os.path.join(r,'_extern.inf'),
                                'r',
                                encoding='latin1')
        extern = inhandler.readlines()
        inhandler.close()

        for line in extern:
            if '\t' in line:
                k = line.strip().split('\t')[0]
                v = line.strip().split('\t')[-1]

                add2dict(ParamDict, k, v.strip())
                all_colnames.append(k)

    except FileNotFoundError as e:
        print('Could not find _extern.inf for file {}:{}'.format(r, e))

    try:
        inhandler = codecs.open(os.path.join(r,'_HEADER.TXT'),
                                'r',
                                encoding='latin1')
        extern = inhandler.readlines()
        inhandler.close()

        for line in extern:
            match= headerPattern.match(line)
            k, v = match.groups()
            add2dict(ParamDict, k, v.strip())
            all_colnames.append(k)

    except FileNotFoundError as e:
        print('Could not find _HEADER.TXT for file {}:{}'.format(r, e))

    data.append(ParamDict)

all_colnames = list(set(all_colnames))

### sort by colnames ###

# Reorder the df
StandardFirstColumns = ['Rawfile',
                        'Acquired Name',
                        'Start Mass',
                        'End Mass',
                        'Cal Date',
                        'Cal Time',

                        'Capillary',
                        'Capillary (kV)',
                        'Desolvation Gas Flow (L/Hr)',
                        'Desolvation Temperature (°C)',
                        'Wave Velocity (m/s)',
                        'Sampling Cone',
                        'Source Temperature (°C)',
                        'Source Wave Height (V)',
                        'Source Wave Velocity (m/s)',
                        'Extraction Cone',
                        'Backing',
                        'Aperture2',

                        'Collision Energy',
                        'Transfer Collision Energy',
                        'Trap Collision Energy',

                        'Trap DC Bias',
                        'Trap DC Entrance',
                        'Trap DC Exit',
                        'Trap Extract Height (V)',
                        'Trap Gas Flow (mL/min)',
                        'Trap Manual Control',
                        'Trap Height (V)',
                        'Trap Wave Height (V)',
                        'Trap Wave Velocity (m/s)',

                        'Transfer DC Entrance',
                        'Transfer DC Exit',
                        'Transfer Extract Height (V)',
                        'Transfer Manual Control',
                        'Transfer Trap Height (V)',
                        'Transfer Wave Height (V)',
                        'Transfer Wave Velocity (m/s)',

                        'Variable Wave Height Enabled',
                        'Variable Wave Velocity Enabled',
                        'Cell pressure [e^-3 mbar]',

                        'IMS',
                        'IMS DC Entrance',
                        'IMS DC Exit',
                        'IMS Extract Height (V)',
                        'IMS Gas Flow (mL/min)',
                        'IMS Manual Control',
                        'IMS Trap Height (V)',
                        'IMS Wave Height (V)',

                        'Detector'
]

fieldnames = list()
for fc in StandardFirstColumns:
    if fc in all_colnames:
        fieldnames.append(fc)

remaining_cols = list()
for c in all_colnames:
    if c not in fieldnames:
        remaining_cols.append(c)

fieldnames.extend(remaining_cols)

d = list()
for idx, c in enumerate(fieldnames):
    col = ['']*len(rawfiles)
    for r, raw in enumerate(rawfiles):
        if c in data[r].keys():
            col[r] = data[r][c]

    d.append(col)

### Optional: Remove identical columns

if args.all is False:
    tmp_header = list()
    tmp = list()
    for idx, col in enumerate(d):
        if len(set(col)) != 1:
            tmp.append(col)
            tmp_header.append(fieldnames[idx])

    d = tmp
    fieldnames = tmp_header

with codecs.open(os.path.join(args.outpath, outname),
                 'w',
                 encoding='utf-8') as out:
    
    writer = csv.writer(out,
                    dialect='excel')

    writer.writerow(fieldnames)
    writer.writerows(list(map(list, zip(*d))))