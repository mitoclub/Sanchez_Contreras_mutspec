#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 12:33:28 2022

@author: scottrk
"""
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
from GlobalVars_ import tissue_type, tissue_type_long, young_mouse_id, color_cycle
from compile_data import calc_depth_plot_data


def get_masked_coords(masked_file):
    with open(masked_file, 'r') as masked:
        masked_coords = []

        for line in masked:
            splitline = line.strip('\n').split('\t')
            masked_coords.append([splitline[1], splitline[2]])

    return masked_coords


def depth_plot(depth_data, uconfint, lconfint, masked_coords, ax, color):
    ax.margins(x=0)
    ax.fill_between(x=range(1, 16300),
                    y1=0,
                    y2=depth_data,
                    color=color
                    )

    ax.fill_between(x=range(1, 16300),
                    y1=uconfint,
                    y2=lconfint,
                    color='black',
                    alpha=0.2
                    )

    for region in masked_coords:
        ax.fill_between(x=range(int(region[0]), int(region[1])), y1=0,
                        y2=28000, alpha=0.2, color='black')

    return ax


def setup_figure():
    fig, axs = plt.subplots(nrows=8, sharey=False, sharex=True,
                            figsize=(8, 10), gridspec_kw=dict(hspace=0.4))

    rcParams['axes.formatter.limits'] = (-6, 6)

    fig.suptitle("Young Mouse Depth", y=0.93, fontsize=12)

    fig.supylabel("Final Duplex Consensus Depth", fontsize=12)

    return fig, axs


if __name__ == "__main__":

    fig, axs = setup_figure()

    masked_coords = get_masked_coords("data/misc_items/mm10_mtDNA_masked_regions.bed")

    for index, tissue in enumerate(tissue_type[:-1]):

        if not os.path.isfile("data/imported_data/Supplemental_Figure_1_" + tissue + "_depth_data.csv"):
            if not os.path.isdir("data/stats/"):
                os.mkdir("data/stats/")

            depth_df = calc_depth_plot_data(young_mouse_id,
                                            tissue,
                                            "data/imported_data/Supplemental_Figure_1_" + tissue + "_depth_data.csv")
        else:
            depth_df = calc_depth_plot_data(young_mouse_id, tissue)

        depth_plot(depth_df['Mean'], depth_df['Std_dev_upper'],
                   depth_df['Std_dev_lower'], masked_coords,
                   axs[index], color_cycle[index])

        axs[index].set_ylim((0, 28000))

        axs[index].set_title(tissue_type_long[index], fontsize=12)

        axs[index].tick_params(labelsize=10)

    axs[7].set_xlabel('Genome Position', size=12)

    if not os.path.isdir("figures"):
        os.mkdir("figures/")

    fig.savefig('figures/Supp_Figure_1.png', dpi=600, facecolor='white', 
                bbox_inches='tight')
