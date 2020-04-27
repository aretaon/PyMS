import pandas as pd
import matplotlib.pylab as plt
import math
import os
from cycler import cycler
infile = r'E:\02_experiments\01_syt1\20171016_jb01_syt1_ddm_membrane_akta\20171018_syt1_sec.asc'

params_to_plot = ['Cond%',
                  #'Conc',
                  ]

data = []
plot_labels = []

with open(infile, 'r') as f:
    f.readline()
    # remove empty entries from the row
    superheader = [x for x in f.readline().strip().split('\t') if x]
    header = f.readline().split('\t')
    header = [x.strip() for x in header]
    for i in range(len(superheader)):
        data.append([[],[]])
    for line in f:
        current = [x.strip() for x in line.split('\t')]
        for idx, el in enumerate(current):
            group = math.floor(idx/ 2) 
            subgroup = idx % 2
            try:
                fl = float(el)
            except:
                if el == '':
                    continue
                else:
                    fl = el
            data[group][subgroup].append(fl)

for h in superheader:
    if 'UV1_280nm' in h:
        uv_280 = h
    elif 'Fractions' in h:
        fractions = h
    else:
        for param in params_to_plot:
            if param in h:
                plot_labels.append(h)


def get_data(label):
    group = superheader.index(label)
    column = data[group]
    return column[0], column[1]

### ´Plot primary y-axis
fig, ax1 = plt.subplots()

lns = []

x_axis, y_axis = get_data(uv_280)

#y_min = min(y_axis)
y_min = -10
y_max = max(y_axis)

x_min = min(x_axis)
# x_min = 40
x_max = max(x_axis)

cur_plot = ax1.plot(x_axis,
                    y_axis,
                    label=uv_280)

lns.extend(cur_plot)

ax1.set_xlabel('Volume [ml]')
ax1.set_ylabel('Absorption at 280 nm')
ax1.set_ylim(0.8 * y_min, 1.2 * y_max)


### plot secondary y-axis
ax2 = ax1.twinx()

ax2.set_prop_cycle(cycler('color', ['c', 'm', 'y', 'k']) +
                   cycler('lw', [1, 2, 3, 4]))

for label in plot_labels:
                  
    x_axis, y_axis = get_data(label)
    y_max = max(y_axis)
    y_min = min(y_axis)
    y_axis = [100 * (y-y_min)/(y_max-y_min) for y in y_axis]
    
    cur_plot = ax2.plot(x_axis,
                        y_axis,
                        label=label)

    lns.extend(cur_plot)

ax2.set_ylabel('Other measures [%]')
ax2.set_ylim(-5, 105)
ax2.set_xlim(x_min, x_max)

### Plot fractions

fr_vol, fr_labels = get_data(fractions)

for v, text in zip(fr_vol, fr_labels):
    if v > x_min and v < x_max:
        ax2.vlines(v, -5, 0, color='r')
        ax2.text(v+1,0, text.strip('"'), color = 'r', rotation = 90)
fig.tight_layout()


### set global things

labs = [l.get_label() for l in lns]
ax2.legend(lns, labs, loc='best')

plt.title('Chromatogram of {}'.format(os.path.basename(infile)))

plt.savefig(os.path.splitext(infile)[0] + '.pdf')