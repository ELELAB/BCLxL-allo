#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.font_manager as fm

import numpy as np
import pandas as pd
import seaborn as sns

def parse_xvg(filename):  
    """Parse the .xvg file derived from the GROMACS g_dist 
    command  and return data as a Pandas DataFrame.
    """
    
    data = []
    with open(filename, "r") as f:
        for line in f:
            # skip description lines
            if (not line.startswith("@") \
                and not line.startswith("#")):          
                
                frame, dist, x, y, z = [float(val.strip("\n")) \
                                        for val in line.split(" ") \
                                        if not val == ""]

                data.append([frame, dist, x, y, z])

    # convert data into a Pandas DataFrame
    data = pd.DataFrame(data, \
                        columns = ["frame", "dist", "x", "y", "z"], \
                        )

    # guess the frame step by subtracting the number of the first
    # frame to the number of the second frame
    frame_step = data["frame"][1] - data["frame"][0]
    
    # add a column storing, for each frame, the percentage of the
    # trajectory covered until that frame
    data["perc"] = \
        ((data["frame"] + frame_step) / data.shape[0])*frame_step

    return data


def plot_gdist(xvg,
               out_name = "g_dist.png", \
               cumulative = True,
               xlabel = u"Distance (Å)", \
               ylabel = "# of frames", \
               xmin = 0.0, \
               xmax = 12.5, \
               xstep = 0.5, \
               ticklabels = True, \
               barcolor = "skyblue", \
               edgecolor = "black", \
               fp_axislabels = None, \
               fp_ticklabels = None, \
                ):
    """Plot data from a .xvg file output by the GROMACS g_dist
    command as a histogram.
    """

    fig, ax = plt.subplots()

    data = parse_xvg(xvg)

    dists = (data["dist"]*10).tolist()

    n, bins, patches = plt.hist(dists, \
                                cumulative = cumulative, \
                                bins = \
                                    np.arange(xmin,xmax+0.5,xstep)-0.25)

    # set the edge color and face color for each bar
    for index, patch in enumerate(patches):
        patch.set_edgecolor(edgecolor)
        patch.set_facecolor(barcolor)
        patch.set_linewidth(0.5)

    # set axis labels
    if fp_axislabels is not None:
        # apply the custom font properties if provided
        ax.set_xlabel(xlabel, \
                      fontproperties = fp_axislabels, \
                      )
        ax.set_ylabel(ylabel, \
                      fontproperties = fp_axislabels, \
                      )
    else:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    
    ax.set_xlim(xmin, xmax)
    ax.set_xticks(np.arange(xmin, xmax+0.5, xstep))

    ax.set_ylim(0,len(dists)+len(dists)/20)
    ax.set_yticks(np.arange(0, len(dists), len(dists)/10))

    # do not show tick labels if 'ticklabels' is False
    if not ticklabels:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    else:
        # generate ticks labels
        yticklabels = ["{:d} %".format(num) \
                       for num in np.arange(0, 101, 10)]

        xticklabels = np.arange(xmin, xmax+0.5, xstep)
        
        if fp_ticklabels is not None:
            # apply the custom font properties if provided
            ax.set_xticklabels(xticklabels, \
                               fontproperties = fp_ticklabels,
                               )
            
            ax.set_yticklabels(yticklabels, \
                               fontproperties = fp_ticklabels,
                               )
        else:
            ax.set_xticklabels(xticklabels)         
            ax.set_yticklabels(yticklabels)

    # save the figure
    plt.savefig(out_name, \
                dpi = 600)


if __name__ == "__main__":

    # relative path to the file containing the g_dist data
    filepath = "../{:s}/g_dist/H73_T75/dist_H73_T75_sidechains.xvg"

    # names of the different ensembles
    names = ["bclxl_puma_3.5nm", "bclxl_puma_capped", \
             "free_2nm_tip3p", "free_2nm_tips3p", "free_cap_tip3p"]

    # get the files containing the g_dist data
    files = [path.format(name) for name in names]

    # get the color palette to be used
    palette = sns.color_palette("colorblind", 5)
    palette2 = palette[0:3] + [palette[4]] + [palette[3]]
    
    # get the path to the font to be used
    font_path = "../fonts/Helvetica-Light.ttf"
    # get the corresponding FontProperties object
    font_prop_labels = fm.FontProperties(fname = font_path, size = 4)

    # for each file
    for xvg, barcolor, name in zip(files, palette2, names):
        
        # set the output name
        out_name = f"gdist_{name}.png"
        
        # generate the plot
        plot_gdist(xvg = xvg, \
                   out_name = out_name, \
                   cumulative = True, \
                   xlabel = "", \
                   ylabel = "", \
                   ticklabels = True, \
                   barcolor = barcolor, \
                   fp_ticklabels = font_prop_labels, \
                   xmax = 20.5)