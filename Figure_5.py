#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 07:59:17 2022

@author: scottrk
"""
import os
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
from GlobalVars_ import mut_type, mut_type_pretty
from compile_data import summary_import, melt_summary


def import_data():
    if not os.path.isfile("data/imported_data/summary_data_tidy.csv"):
        if not os.path.isdir("data/imported_data/"):
            os.mkdir("data/imported_data")

        data = melt_summary(summary_import("data/Mouse_aging_mtDNA_summary.csv"))
        data.to_csv("data/imported_data/summary_data_tidy.csv")
    else:
        data = pd.read_csv("data/imported_data/summary_data_tidy.csv")
    
    data = data.query("Treatment!='perf'& Age=='Old' & Tissue in ['K', 'M', 'He'] & \
                       Class not in ['Total_SNV_Freq', 'Total_InDel_Freq']")
    
    return data


def setup_figure():
    
    mosaic = '''AABB
                .CC.'''
                
    fig = plt.figure(constrained_layout=True, figsize=(12, 8))
    ax = fig.subplot_mosaic(mosaic, gridspec_kw=dict(hspace=0.1))
    subax = inset_axes(ax['A'], width="17%", height="40%", 
                       bbox_to_anchor=[-0.31, -0.15, 1, 1], bbox_transform=ax['A'].transAxes)
    
    subax.margins(y=0)
    ax['A'].set_ylim(0, 1.6e-5)
    ax['A'].set_title("Kidney", fontsize=20)
    ax['B'].set_ylim(0, 2.8e-6)
    ax['B'].set_title("Heart", fontsize=20)
    ax['C'].set_ylim(0, 4e-6)
    ax['C'].set_title("Sk. Muscle", fontsize=20)
    
    return fig, ax, subax


if __name__ == '__main__':    
    
    data = import_data()
    fig, ax, subax = setup_figure()
    subplot_id = {0: 'A', 1: 'B', 2: 'C'}
    
    for i, tissue in enumerate(['K', 'He', 'M']):
        
        j = subplot_id[i]
        
        sns.barplot(x="Class", y="Frequency", hue="Treatment", 
                    data=data.query("Tissue==@tissue"), order=mut_type, ci='sd', 
                    edgecolor='black', lw=1.5, errwidth=1.5, capsize=0.1, 
                    errcolor='black', ax=ax[j])
        sns.despine(ax=ax[j])
        
        ax[j].margins(y=0)
        fig.canvas.draw()
        plt.setp(ax[j].get_yaxis().get_offset_text(), visible=False)
        order_of_mag = ax[j].get_yaxis().get_offset_text().get_text()[-2:]
        string = "Mutation Frequency ($\mathregular{10^{" + str(order_of_mag) + "}}$)"
        ax[j].set_xticklabels(mut_type_pretty, rotation=45)
        ax[j].set_ylabel(string, fontsize=16)
        ax[j].set_xlabel("")
        ax[j].tick_params(labelsize=14)
        ax[j].legend(fontsize=14)
        
        hatches = ['///',  '---']
        for x, bar in enumerate(ax[j].patches):
            ax[j].patches[x].set_facecolor(('#ff9200', '#ff2600', '#008e00')[i])
            if x in range(6, 12):
                bar.set_hatch(hatches[0])
            elif x in range(12, 18):
                bar.set_hatch(hatches[1])
        ax[j].legend(markerscale=5, fontsize=16)        

    sub_data = data.query("Tissue=='K' & Class=='C>G/G>C' & Treatment!='perf'")   
    
    sns.barplot(x='Class', y='Frequency', hue='Treatment', data=sub_data,
                palette='bright', ci='sd', edgecolor='black', lw=1.5,
                errwidth=1.5, capsize=0.04, errcolor='black', ax=subax)
    
    hatches = ['///',  '---']
    for x, bar in enumerate(subax.patches):
        subax.patches[x].set_facecolor('#ff9200')
        if x == 1:
            bar.set_hatch(hatches[0])
        elif x == 2:
            bar.set_hatch(hatches[1])
    
    fig.canvas.draw()
    plt.setp(subax.get_yaxis().get_offset_text(), visible=False)
    order_of_mag = subax.get_yaxis().get_offset_text().get_text()[-2:]
    string = "Mut. Freq ($\mathregular{10^{" + str(order_of_mag) + "}}$)"
    subax.set_ylabel(string, fontsize=15)
    x = subax.get_xticklabels()
    subax.set_xticklabels(['G→C/C→G'])
    subax.set_xlabel('')
    subax.legend_.remove()
    subax.set_ylim(0, 1.3e-6)
    subax.tick_params(labelsize=14)
    
    if not os.path.isdir('figures/'):
        os.mkdir('figures/')
        
    fig.savefig("figures/Figure_5.png", facecolor='white', dpi=600) 
    fig.savefig("figures/Figure_5.pdf", facecolor='white', dpi=600) 
