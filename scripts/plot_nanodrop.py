# -*- coding: utf-8 -*-
"""
Created on Fri May  8 18:40:40 2020

@author: User
"""
import matplotlib.pyplot as plt
import datetime
import argparse

parser = argparse.ArgumentParser(description='Generate Figure from nanodrop UV-Vis data')
parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
parser.add_argument('--save', default=False)

args = parser.parse_args()

data = []

for f in args.file:
    header = [x.strip('"\n') for x in f.readline().split(',')]
    for line in f.readlines():
        l = [x.strip('"\n') for x in line.split(',')]
        data.append(l)

fig = plt.figure()
ax = fig.gca()

uvrange = [float(x) for x in header[20:]]

for el in data:
    if el[1] == 'Blank':
        continue
    signal = [float(x) for x in el[20:]]
    ax.plot(uvrange, signal, label="{} {} mg/mL".format(el[1], el[2]))

plt.xlabel('Wavelength [nm]')
plt.ylabel('Absorbance')

plt.legend()

if args.save:
    plt.savefig('{}_uv_absorbance.pdf'.format(str(datetime.datetime.now().date())))
else:
    plt.show()