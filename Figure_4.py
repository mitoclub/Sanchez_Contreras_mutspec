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
from GlobalVars_ import color_cycle, tissue_type, tissue_type_long, tissue_type_abbrev, \
                        mut_type_pretty
from compile_data import mut_file_import, calc_clone_numbers


def setup_figure():
    rcParams['axes.formatter.limits'] = (-6, 6)
    fig, ax = plt.subplots(ncols=2, figsize=(25, 8))

    ax[0].axhline(0.0005, color='black', ls='--', lw=2, alpha=1)
    ax[1].axhline(0.0005, color='black', ls='--', lw=2, alpha=1)

    ax1ins = inset_axes(ax[1], width="60%", height="60%",
                        bbox_to_anchor=[-0.05, -0.05, 1, 1],
                        bbox_transform=ax[1].transAxes)

    ax1ins.set_ylim(0, 0.0005)

    return fig, ax, ax1ins

def spectrum(data, mut_type, ax, fill=False, ymin=None, ymax=None, legend=True):
    
    sns.barplot(x='Class', y='Frequency', hue='Tissue', data=data,
                       order=mut_type, ci='sd',
                       edgecolor='black', lw=1.5, errwidth=1.7,
                       capsize=0.07, errcolor='black', ax=ax)

    sns.stripplot(x="Class", y="Frequency", hue="Tissue", data=data,
                  order=mut_type, ax=ax, alpha=0.7, dodge=True, color='black')

    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

    if ymin is not None and ymax is not None:
        ax.set_ylim(ymin, ymax)
    
    patch_list = []
    
    if not fill:
        for i, _ in enumerate(ax.patches):

            j = i // 6
            ax.patches[i].set_facecolor(color_cycle[j])
            r, g, b, a = ax.patches[i].get_facecolor()
            ax.patches[i].set_facecolor((r, g, b, 0.15))  
            ax.patches[i].set_edgecolor(color_cycle[j])
            ax.patches[i].set(lw=2.4)

            if i % 6 == 0:
                patch_list.append([r, g, b, 0.15])
        
        if legend:
            legend = [Patch.Patch(facecolor=patch_list[x], lw=2.4, edgecolor=color_cycle[x],
                                  label=tissue_type_long[x]) for x in range(8)]
                
    
    else:

        for i, _ in enumerate(ax.patches):

            j = i // 6
            ax.patches[i].set_facecolor(color_cycle[j])
            ax.patches[i].set_edgecolor('black')

            if legend:
                legend = [Patch.Patch(facecolor=color_cycle[x], edgecolor='black',
                                      label=tissue_type_long[x]) for x in range(8)]

    ax.set_xticklabels(mut_type_pretty)
    ax.set_xlabel('')
    ax.tick_params(labelsize=18)
    ax.margins(x=.01)

    ax.set_ylabel("Clone Frequency", fontsize=25)
    ax.set_xlabel("")
    ax.set_xticklabels(mut_type_pretty, rotation=45)
    ax.tick_params(labelsize=22)

    return ax, legend

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

    plot_A, legend0 = spectrum(data=clone_data_long.query("Cohort=='Young'"), 
                               mut_type=mut_type_names, ax=ax[0], fill=False)

    plot_B, legend1 = spectrum(data=clone_data_long.query("Cohort=='Old'"), 
                               mut_type=mut_type_names, ax=ax[1], fill=True)
    
    sns.despine(ax=ax[0])
    sns.despine(ax=ax[1])

    plot_B_ins, _ = spectrum(data=clone_data_long.query("Cohort=='Old'"), 
                                   mut_type=mut_type_names, ax=ax1ins, legend=False,
                                   fill=True, ymin=0, ymax=0.0005)    


    ax[0].legend(handles=legend0, fontsize=22, ncol=4, bbox_to_anchor=[1.1, 1.2], frameon=False)
    ax[1].legend(handles=legend1, fontsize=22, ncol=4, bbox_to_anchor=[1.2, 1.2], frameon=False)



    ax[0].set_ylim(0, 0.0011)
    ax[1].set_ylim(0, 0.0125)
    

    ax1ins.set_xticklabels(mut_type_pretty, fontsize=17, rotation=45)

    ax1ins.set_xlabel("")

    ax1ins.set_ylabel("Clone Frequency", fontsize=18)
    ax1ins.tick_params(labelsize=16)
    ax1ins.legend_.remove()

    if not os.path.isdir("figures"):
        os.mkdir("figures/")

    fig.savefig("figures/Figure_4A-B.png", facecolor='white', dpi=600, 
                bbox_inches='tight')
    fig.savefig("figures/Figure_4A-B.pdf", facecolor='white', dpi=600, 
                bbox_inches='tight')
