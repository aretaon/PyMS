import matplotlib.pyplot as plt
import matplotlib.font_manager
import matplotlib.patches

import sys

## Text and alignment

def PlotXlinkMap(pepseq1,
                 pepseq2,
                 xlink1,
                 xlink2,
                 prec_ch,
                 positions,
                 types,
                 charges,
                 peptide_types,
                 fontsize=12,
                 max_aa_per_window = 20,
                 label_origin = [0,0],
                 ax=None):
                   
    global fp
    
    fp = matplotlib.font_manager.FontProperties(
                    family='monospace', style='normal', size=12,
                    weight='normal', stretch='normal')

    if ax is None:
        ax = plt.gca()
        
    fig = plt.gcf()
    r = fig.canvas.get_renderer()
        
    # set the limits of the figure
    ax.set_ylim([-0.5, 0.5])
    ax.set_xlim([-0.5,0.5])
    
    yoffset = 0.2
    
    xoffset = abs(xlink1-xlink2) / 2
    
    def draw_text(pepseq, yoffset, xoffset=0, verticalalignment='center', ax=None):
    
        peporigin = -0.5*len(pepseq1)/max_aa_per_window + (xoffset)/max_aa_per_window
        x = peporigin 
        xstep = 1/max_aa_per_window
    
        if ax is None:
            ax = plt.gca()
    
        for idx, aa in enumerate(pepseq):
            y = yoffset
            s = aa
    
            t = ax.text(x, y, s,
                        verticalalignment=verticalalignment,
                        fontproperties = fp
                        )
            print(x,y,s)
            x += xstep
            
        return peporigin
    
    ori1 = draw_text(pepseq1, yoffset, -xoffset, ax=ax)
    ori2 = draw_text(pepseq2, -yoffset, xoffset, ax=ax)
    
    x1 = ori1 + (xlink1-0.75)/max_aa_per_window
    x2 = ori2 + (xlink2-0.75)/max_aa_per_window
    ax.plot((x1, x2), (-yoffset*0.8, yoffset*0.8), color='lightgrey', lw=3)
 
    ## Annotation
    
    oris = (ori1, ori2)
    
    def annotate_fragment(aa_lost, fragment_type, charge, peptide='alpha', oris=oris, ax=ax, yoffset=yoffset):
        if peptide == 'alpha':
            ori = oris[0]
            label_sign = 1
        elif peptide == 'beta':
            ori = oris[1]
            label_sign = -1
            yoffset = -yoffset
        else:
            return
        
        nterm = False
        cterm = False
        
        if fragment_type in ['a', 'A', 'b', 'B', 'c', 'C']:
            nterm = True
        elif fragment_type in ['x', 'X', 'y', 'Y', 'z', 'Z']:
            cterm = True
        
        x = ori + (aa_lost-0.33)/max_aa_per_window
        ax.plot((x, x), (0.8*yoffset, 1.2*yoffset), 'k-')
    
        # charges
        if charge > 4:
            charge = 4
        
        cmap_nterm = ['lightsalmon', 'tomato', 'r', 'brown']
        cmap_cterm = ['lightblue', 'deepskyblue', 'darkturquoise', 'darkcyan']
        
        if nterm:
            label_color = cmap_nterm[charge-1]
        elif cterm:
            label_color = cmap_cterm[charge-1]
    
        # arrows for fragment direction
        if nterm:
            ax.plot((x, x-0.25/max_aa_per_window), (1.2*yoffset, 1.3*yoffset), 'k-')
        elif cterm:
            ax.plot((x, x+0.25/max_aa_per_window), (0.8*yoffset, 0.7*yoffset), 'k-')
    
        # labels for the charges
        if nterm:
            yposition =  1.4*yoffset + label_sign * (charge-1)*0.03
            ax.plot(x, yposition, marker='o', markersize=4, color=label_color)
        elif cterm:
            yposition =  0.6*yoffset - label_sign * (charge-1)*0.03
            ax.plot(x, yposition, marker='o', markersize=4,color=label_color)

    for idx, pos in enumerate(positions):
        annotate_fragment(pos, types[idx], charges[idx], peptide=peptide_types[idx])

    ax.axis('off')