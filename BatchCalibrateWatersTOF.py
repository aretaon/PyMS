# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 13:43:53 2019

@author: User
"""

def getKeyAndValue(line):
    key = line.split('=')[0].strip()
    value = line.split('=')[1].strip()
    return key, value

sclEntries = dict()

with open(r'G:\01_data\Synapt\20190116_ADH_CONA_IMS\20190116_JB_CSI_cal.scl', 'r') as f:
    for line in f.readlines():
        # skip heading lines
        if not line.startswith('['):
            key, value = getKeyAndValue(line)
            sclEntries[key] = value

