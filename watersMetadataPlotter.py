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

args = parser.parse_args()

#args.dir = r'C:\Users\User\Documents\03_software\python\PyMS\watersMetadataPlotter_testdata'

rawfiles = list()

for root, dirs, files in os.walk(args.dir):
    for f in dirs:
        if f.endswith('.raw'):
            rawfiles.append(os.path.join(root, f))

base_rawfiles = [os.path.split(f)[1][:-4] for f in rawfiles]

outname = '{:s}_complete_tune.csv'.format(os.path.basename(args.dir))

DictList = []

headerPattern = re.compile(r'\$\$ ([^:]+)\:(.*)')

for r in rawfiles:

    ParamDict = {}

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
                key = line.strip().split('\t')[0]
                value = line.strip().split('\t')[-1]
                ParamDict[key] = value

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
            key, value = match.groups()
            ParamDict[key.strip()] = value.strip()

    except FileNotFoundError as e:
        print('Could not find _HEADER.TXT for file {}:{}'.format(r, e))


    DictList.append(ParamDict)

fieldnames = set()

for d in DictList:
    # use the union parameter to add the dict to the fieldnames set
    fieldnames |= set(d.keys())

# generate an order to print
fieldnames = sorted(list(fieldnames))


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


def MoveForward(name, position, data):
    """
    Move a list entry at a specified position
    
    Requires:
    - name: Name of entry header
    - position: Desired final position
    - data: Pandas dataframe
    
    Returns:
    - ordered dataframe
    """
    
    rawfile_idx = data.index(name)
    # move the rawfile ID at the indicated position
    data[position], data[rawfile_idx] = \
    data[rawfile_idx], data[position]

    return data


outhandler = codecs.open(os.path.join(args.outpath, outname),
                         'w',
                         encoding='utf-8')

firstEntries = list()

for name in StandardFirstColumns:
    if name in fieldnames:
        firstEntries.append(name)

for name, position in zip(firstEntries, range(len(firstEntries))):
    # names differ depending on the MS-Instrument used
    fieldnames = MoveForward(name, position, fieldnames)

writer = csv.DictWriter(outhandler,
                        fieldnames=fieldnames,
                        dialect='excel')
writer.writeheader()
for d in DictList:
    writer.writerow(d)
    
outhandler.close()