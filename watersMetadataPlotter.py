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
import time
import codecs
import csv
import re

thisdate = time.strftime("%Y%m%d")

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--dir", default=os.getcwd(), help="Path to the directory containing the raw files")

parser.add_argument("-o", "--outpath", default=os.getcwd(), help="Path to a folder to write the output")

parser.add_argument('-a', "--all", action="store_true", default=False, help="Print all columns and not only those that differ between the samples")

args = parser.parse_args()

#args.dir = r'C:\Users\User\Documents\03_software\python\PyMS\watersMetadataPlotter_testdata'
#args.all = False

rawfiles = list()

for root, dirs, files in os.walk(args.dir):
    for f in dirs:
        if f.endswith('.raw'):
            rawfiles.append(os.path.join(root, f))

base_rawfiles = [os.path.split(f)[1][:-4] for f in rawfiles]

outname = '{:s}_complete_tune.csv'.format(os.path.basename(args.dir))

colnames2values = dict()

headerPattern = re.compile(r'\$\$ ([^:]+)\:(.*)')

for r in rawfiles:

    ParamDict = {}
    
    keys = list()
    values = list()

    ### Collection data ###
    ParamDict['Rawfile'] = os.path.split(r)[1]
    # encoding has to be set to latin1 as MassLynx does not use UTF-8
    try:
        inhandler = codecs.open(os.path.join(r,'_extern.inf'),
                                'r',
                                encoding='latin1')
        extern = inhandler.readlines()
        inhandler.close()

        for line in extern:
            if '\t' in line:
                keys.append(line.strip().split('\t')[0])
                values.append(line.strip().split('\t')[-1])
#                ParamDict[key] = value

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
            keys.append(k.strip())
            values.append(v.strip())
#            ParamDict[key.strip()] = value.strip()

    except FileNotFoundError as e:
        print('Could not find _HEADER.TXT for file {}:{}'.format(r, e))

    ### Appending data to datatable

    for idx, k in enumerate(keys):
        if k in colnames2values.keys():
            colnames2values[k].append(values[idx])
        else:
            colnames2values[k] = [values[idx]]


### Optional: Remove identical columns

if args.all is False:
    noninformative_cols  =list()
    for colname in colnames2values.keys():
        if len(set(colnames2values[colname])) == 1:
            noninformative_cols.append(colname)

    for k in noninformative_cols:
        del colnames2values[k]
        

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
    if fc in colnames2values.keys():
        fieldnames.append(fc)

remaining_cols = list()
for c in colnames2values.keys():
    if c not in fieldnames:
        remaining_cols.append(c)
        
fieldnames.extend(remaining_cols)

outhandler = codecs.open(os.path.join(args.outpath, outname),
                         'w',
                         encoding='utf-8')

writer = csv.writer(outhandler,
                    dialect='excel')

writer.writerow(fieldnames)
writer.writerows(zip(*[colnames2values[key] for key in fieldnames]))

outhandler.close()