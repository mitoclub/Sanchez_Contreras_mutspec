#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:00:19 2022

@author: scottrk
"""
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from GlobalVars_ import tissue_type, tissue_type_long, color_cycle, mut_type, \
     mut_type_pretty, tissue_type_abbrev
from compile_data import summary_import, melt_summary


def import_data():

    if not os.path.isfile("data/imported_data/summary_data_tidy.csv"):
        if not os.path.isdir("data/imported_data/"):
            os.mkdir("data/imported_data")

        data = melt_summary(summary_import("data/Mouse_aging_mtDNA_summary.csv"))
        data.to_csv("data/imported_data/summary_data_tidy.csv")
    else:
        data = pd.read_csv("data/imported_data/summary_data_tidy.csv")

    data = data.query("Treatment!='perf'& Age=='Old' \
                      & Class not in ['Total_SNV_Freq', 'Total_InDel_Freq']")

    return data


def setup_figure():

    mosaic = '''AABB
                CCDD'''

    fig = plt.figure(constrained_layout=True, figsize=(12, 8))
    ax = fig.subplot_mosaic(mosaic, gridspec_kw=dict(hspace=0.1, wspace=0.2))
    ax['A'].set_ylim(0, 1e-5)
    ax['B'].set_ylim(0, 1e-5)
    ax['C'].set_ylim(0, 8e-6)
    ax['D'].set_ylim(0, 5e-6)

    return fig, ax


if __name__ == '__main__':

    data = import_data()
    fig, ax = setup_figure()
    subplot_id = {0: 'A', 1: 'B', 2: 'C', 3: 'D'}
    intervention_color = color_cycle[2:6]

    for i, tissue in enumerate(tissue_type_abbrev[2:6]):

        j = subplot_id[i]

        plot = sns.barplot(x="Class", y="Frequency", hue="Treatment",
                           data=data.query("Tissue==@tissue"), order=mut_type,
                           ci='sd', edgecolor='white', lw=2, errwidth=1.5,
                           capsize=0.1, errcolor='black', ax=ax[j])

        ax[j].margins(y=0)
        fig.canvas.draw()
        plt.setp(ax[j].get_yaxis().get_offset_text(), visible=False)
        order_of_mag = ax[j].get_yaxis().get_offset_text().get_text()[-2:]
        string = "Mutation Frequency ($\mathregular{10^{" + str(order_of_mag) + "}}$)"
        ax[j].set_ylabel(string, fontsize=16)
        ax[j].set_xlabel('')
        ax[j].set_xticklabels(mut_type_pretty, rotation=45)
        ax[j].tick_params(labelsize=14)
        ax[j].set_title(tissue_type_long[2:6][i], fontsize=16, y=0.85)

        hatches = ['///', '--']
        for x, bar in enumerate(ax[j].patches):
            plot.patches[x].set_facecolor(intervention_color[i])
            if x in range(6, 12):
                bar.set_hatch(hatches[0])
            elif x in range(12, 18):
                bar.set_hatch(hatches[1])
        plot.legend(markerscale=5, fontsize=16)

    if not os.path.isdir('figures/'):
        os.mkdir('figures/')

    fig.savefig("figures/Supp_Figure_12.png", facecolor='white', dpi=600)
    fig.savefig("figures/Supp_Figure_12.pdf", facecolor='white', dpi=600)
