#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 13:07:43 2022

@author: scottrk
"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as Patch
from GlobalVars_ import color_cycle, tissue_type, tissue_type_long, tissue_type_abbrev
from HelperFuncs_ import p_val_convert
from compile_data import melt_summary, summary_import
import os
from compute_stats import heatmap_stats


def setup_figure():
    mosaic = """AA
                BC"""
    mosaic2 = """AA
                 BC"""
    fig = plt.figure(constrained_layout=True, figsize=(15, 8))
    left, right = fig.subfigures(nrows=1, ncols=2, wspace=0.08)
    axd = left.subplot_mosaic(mosaic, gridspec_kw=dict(hspace=0.1))
    axd2 = right.subplot_mosaic(mosaic2, gridspec_kw=dict(hspace=0.1))

    return fig, axd, axd2


def mutation_bar(y, data, ax):
    fc = ['white'] * 8 + color_cycle
    ec = color_cycle + ['black'] * 8

    plot = sns.barplot(x="Tissue", y=y, hue="Age", data=data, order=tissue_type_abbrev[: -1],
                       palette='bright', ci='sd', edgecolor='black', lw=1.2,
                       errwidth=1.5, capsize=0.1, errcolor='black', ax=ax)

    sns.stripplot(x="Tissue", y=y, hue="Age", data=data, order=tissue_type_abbrev[:-1],
                  dodge=True, ax=plot, alpha=0.7, color='black')

    sns.despine(ax=plot)

    for i in range(16):

        ax.patches[i].set_facecolor(fc[i])
        ax.patches[i].set_edgecolor(ec[i])

        if i < 8:
            ax.patches[i].set(lw=2.4)

    plt.setp(ax.get_yaxis().get_offset_text(), visible=False)
    ax.set_ylabel('Mutation Frequency($\mathregular{10^{-6}}$)', fontsize='xx-large')
    ax.set_xlabel('', fontsize='xx-large')
    ax.tick_params('y', labelsize='xx-large')
    ax.tick_params('x', labelsize='x-large')
    ax.set_xticklabels(tissue_type_long, rotation=45, fontdict={'horizontalalignment': 'center'})

    legend = [Patch.Patch(facecolor='white', edgecolor='black', label='Young'),
              Patch.Patch(facecolor='lightgrey', edgecolor='black', label='Old')]

    ax.legend(handles=legend, fontsize='x-large')

    return plot


def mod_heatmap(data, labels, ax):
    plot = sns.heatmap(data, square=True, cbar=False, cmap='rocket_r',
                       linewidths=.5, linecolor='black', vmin=0, vmax=10,
                       xticklabels=labels, yticklabels=labels, ax=ax)

    sns.despine(left=True, bottom=True, ax=plot)
    plot.set_facecolor('silver')
    plot.tick_params(labelsize=18)

    return plot


if __name__ == '__main__':

    if not os.path.isfile("data/imported_data/summary_data_tidy.csv"):
        if not os.path.isdir("data/imported_data/"):
            os.mkdir("data/imported_data")

        data = melt_summary(summary_import("data/Mouse_aging_mtDNA_summary.csv"))
        data.to_csv("data/imported_data/summary_data_tidy.csv")
    else:
        data = pd.read_csv("data/imported_data/summary_data_tidy.csv")

    fig, axd, axd2 = setup_figure()

    mutation_bar("Frequency", data.query("Treatment=='NT' & Class=='Total_SNV_Freq'"), axd['A'])
    mutation_bar("Frequency", data.query("Treatment=='NT' & Class=='Total_InDel_Freq'"), axd2['A'])

    if not os.path.isfile("data/stats/Figure_1A_Young_stats.csv"):
        if not os.path.isdir("data/stats/"):
            os.mkdir("data/stats/")

        young_SNV = heatmap_stats(data, 'Young', 'Total_SNV_Freq', "data/stats/Figure_1A_Young_stats.csv")
    else:
        young_SNV = pd.read_csv("data/stats/Figure_1A_Young_stats.csv",
                                index_col=0)

    if not os.path.isfile("data/stats/Figure_1A_Old_stats.csv"):
        if not os.path.isdir("data/stats/"):
            os.mkdir("data/stats/")

        old_SNV = heatmap_stats(data, 'Old', 'Total_SNV_Freq')
    else:
        old_SNV = pd.read_csv("data/stats/Figure_1A_Old_stats.csv",
                              index_col=0)

    if not os.path.isfile("data/stats/Figure_1D_Young_stats.csv"):
        if not os.path.isdir("data/stats/"):
            os.mkdir("data/stats/")

        young_InDel = heatmap_stats(data, 'Young', 'Total_InDel_Freq')
    else:
        young_InDel = pd.read_csv("data/stats/Figure_1D_Young_stats.csv",
                                  index_col=0)

    if not os.path.isfile("data/stats/Figure_1D_Old_stats.csv"):
        if not os.path.isdir("data/stats/"):
            os.mkdir("data/stats/")

        old_InDel = heatmap_stats(data, 'Old', 'Total_InDel_Freq')
    else:
        old_InDel = pd.read_csv("data/stats/Figure_1D_Old_stats.csv",
                                index_col=0)

    mod_heatmap(p_val_convert(young_SNV), tissue_type_abbrev[:-1], axd['B'])
    mod_heatmap(p_val_convert(old_SNV), tissue_type_abbrev[:-1], axd['C'])
    mod_heatmap(p_val_convert(young_InDel), tissue_type_abbrev[:-1], axd2['B'])
    mod_heatmap(p_val_convert(old_InDel), tissue_type_abbrev[:-1], axd2['C'])

    for i, j in enumerate(['B', 'C']):
        title = ['Young', 'Old']
        axd[j].set_title(title[i], fontsize='xx-large')
        axd2[j].set_title(title[i], fontsize='xx-large')

    if not os.path.isdir("figures/"):
        os.mkdir("figures")

    fig.savefig('figures/Figure_1.png', dpi=600, facecolor='white')
    fig.savefig('figures/Figure_1.pdf', dpi=600, facecolor='white')