#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 09:20:55 2022

@author: scottrk
"""
import os
import matplotlib.pyplot as plt
import matplotlib.patches as Patch
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import pandas as pd
import itertools
from GlobalVars_ import color_cycle, tissue_type, tissue_type_long, mut_type_pretty
from compile_data import mut_file_import, calc_clone_numbers


def setup_figure():
    rcParams['axes.formatter.limits'] = (-6, 6)
    fig, ax = plt.subplots(ncols=2, figsize=(25, 8))

    ax[0].axhline(0.0005, color='black', ls='--', alpha=0.5)
    ax[1].axhline(0.0005, color='black', ls='--', alpha=0.5)

    ax1ins = inset_axes(ax[1], width="60%", height="60%",
                        bbox_to_anchor=[-0.05, -0.05, 1, 1],
                        bbox_transform=ax[1].transAxes)

    ax1ins.set_ylim(0, 0.0005)

    return fig, ax, ax1ins


if __name__ == '__main__':

    mut_type_names = ["G>A/C>T_Freq", "A>G/T>C_Freq", "G>T/C>A_Freq", "G>C/C>G_Freq", "A>T/T>A_Freq", "A>C/T>G_Freq"]

    if not os.path.isfile("data/imported_data/summary_clone_data.csv"):
        if not os.path.isdir("data/imported_data"):
            os.mkdir("data/imported_data/")

        mut_data = mut_file_import()
        final_clone_data = calc_clone_numbers(mut_data)

        mut_data.to_csv("data/imported_data/mut_file_data.csv")
        final_clone_data.to_csv("data/imported_data/summary_clone_data.csv")
    else:
        mut_data = pd.read_csv("data/imported_data/mut_file_data.csv",
                               index_col=[0, 1])
        final_clone_data = pd.read_csv("data/imported_data/summary_clone_data.csv")

    clone_data_long = pd.melt(final_clone_data, id_vars=['Mouse_ID', 'Tissue', 'Cohort'],
                              value_vars=['A>T/T>A_Freq', 'A>C/T>G_Freq', 'A>G/T>C_Freq', 'G>T/C>A_Freq',
                                          'G>A/C>T_Freq', 'G>C/C>G_Freq'])
    clone_data_long.columns = ["MouseID", "Tissue", 'Cohort', 'Class', 'Frequency']

    fig, ax, ax1ins = setup_figure()

    sns.barplot(x='Class', y='Frequency', hue='Tissue',
                data=clone_data_long.query("Cohort=='Young'"),
                hue_order=tissue_type[:-1],
                palette=color_cycle, facecolor='white', order=mut_type_names,
                ci='sd', lw=1.5, errwidth=1.2, errcolor='black',
                capsize=0.07, ax=ax[0])

    sns.despine(ax=ax[0])

    sns.barplot(x='Class', y='Frequency', hue='Tissue',
                data=clone_data_long.query("Cohort=='Old'"),
                hue_order=tissue_type[:-1], palette=color_cycle,
                order=mut_type_names, ci='sd', edgecolor='black', lw=1.2,
                errwidth=1.2, errcolor='black', capsize=0.07, ax=ax[1])

    sns.despine(ax=ax[1])

    sns.barplot(x='Class', y='Frequency', hue='Tissue',
                data=clone_data_long.query("Cohort=='Old'"),
                hue_order=tissue_type[:-1], palette=color_cycle,
                order=mut_type_names, edgecolor='black', ci='sd', lw=1,
                errwidth=1.2, errcolor='black', capsize=0.07, ax=ax1ins)

    ec = zip(color_cycle, color_cycle, color_cycle, color_cycle, color_cycle, color_cycle)
    ec = list(itertools.chain.from_iterable(ec))

    for i in range(len(ax[0].patches)):
        ax[0].patches[i].set_edgecolor(ec[i])
        legend0 = [Patch.Patch(facecolor='white', edgecolor=color_cycle[i],
                               label=tissue_type_long[i]) for i in range(8)]

        legend1 = [Patch.Patch(facecolor=color_cycle[i], edgecolor='black',
                               label=tissue_type_long[i]) for i in range(8)]

    ax[0].legend(handles=legend0, fontsize=14, ncol=4, bbox_to_anchor=[0.9, 1.12], frameon=False)
    ax[1].legend(handles=legend1, fontsize=14, ncol=4, bbox_to_anchor=[0.9, 1.12], frameon=False)

    ax[0].set_ylabel("Clone Frequency", fontsize=18)
    ax[1].set_ylabel("Clone Frequency", fontsize=18)

    ax[0].set_xlabel("")
    ax[1].set_xlabel("")

    ax[0].set_xticklabels(mut_type_pretty, fontsize=18)
    ax[1].set_xticklabels(mut_type_pretty, fontsize=18)

    ax[0].tick_params(labelsize=15)
    ax[1].tick_params(labelsize=15)

    ax[0].set_ylim(0, 0.0011)
    ax[1].set_ylim(0, 0.0125)

    ax1ins.set_xticklabels(mut_type_pretty, fontsize=12, rotation=45)

    ax1ins.set_xlabel("")

    ax1ins.set_ylabel("Clone Frequency", fontsize=14)

    ax1ins.legend_.remove()

    if not os.path.isdir("figures"):
        os.mkdir("figures/")

    fig.savefig("figures/Figure_4.png", facecolor='white', dpi=600, 
                bbox_inches='tight')
