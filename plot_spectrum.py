# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

#header = 'IGNSample Number,IGNSample Name,IGNmg/mL,IGNUnits,IGNSample Type,IGNE1%,IGNA280,IGNPathlength (mm),IGN260/280,IGN260/280 Alert,IGNDate/Time,IGNMW,IGNExt Coeff,IGNBaseline nm,IGNBaseline abs,IGNAccount,IGNSession,IGNMeasure Type,IGNInstrument,IGNExposure,"220","221","222","223","224","225","226","227","228","229","230","231","232","233","234","235","236","237","238","239","240","241","242","243","244","245","246","247","248","249","250","251","252","253","254","255","256","257","258","259","260","261","262","263","264","265","266","267","268","269","270","271","272","273","274","275","276","277","278","279","280","281","282","283","284","285","286","287","288","289","290","291","292","293","294","295","296","297","298","299","300","301","302","303","304","305","306","307","308","309","310","311","312","313","314","315","316","317","318","319","320","321","322","323","324","325","326","327","328","329","330","331","332","333","334","335","336","337","338","339","340","341","342","343","344","345","346","347","348","349","350"'.split(',')
#header = [x.strip('"') for x in header]

infile = r'C:\Users\User\Downloads\20191220_gab1_auswertung\20191213_gab1.csv'

data = pd.read_csv(infile)

data_cols = []
for col in data.columns:
    try:
        data_cols.append(str(int(col)))
    except:
        print('Column {} is not a data-column'.format(col))

data_cols.append('Sample Name')

spectrum = data[data_cols].set_index('Sample Name').transpose()

ax = spectrum.plot()

ax.set_ylabel('Int')
ax.set_xlabel('Wavelength (nm)')

fig = ax.get_figure()

#plt.savefig(os.path.splitext(infile)[0] + '.pdf')