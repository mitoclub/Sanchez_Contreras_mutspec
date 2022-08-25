#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 14:23:13 2022

@author: scottrk
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from GlobalVars_ import tissue_type, tissue_type_long, mut_type, \
                        mut_type_pretty, young_mouse_id, old_mouse_id
                        

def process_data(data, age):
    
    summed_data = []
    summed_data = []
    #df = pd.DataFrame()

    for tissue in tissue_type[:-1]:
        
        count_cols = [x + "_Count" for x in mut_type]
        columns = count_cols + ['Total_SNVs']
        sub_data = data.query("Tissue==@tissue & Treatment=='NT' & Age==@age")[columns].sum(axis=0)
        summed_data.append(pd.DataFrame(sub_data))

    df = pd.concat(summed_data, axis=1)
    df.columns = tissue_type[:-1]
    df = df.T
    df[mut_type_pretty] = df[[x + "_Count" for x in mut_type]].div(df['Total_SNVs'].mul(0.01), axis=0)
    
    return df 


def plot_figure():
    
    fig, ax = plt.subplots(nrows=2, figsize=(12,8), gridspec_kw={'hspace':0.8})
    colors = ['#fffb00', '#d4fb79', '#ff40ff', '#ff85ff', '#73fa79', '#019051']
    cmap = LinearSegmentedColormap.from_list("mycmap", colors)
    
    young[mut_type_pretty].plot(kind='bar', stacked=True, ec='black', lw=0.4,
                                ax=ax[0], colormap=cmap)
    old[mut_type_pretty].plot(kind='bar', stacked=True, ec='black', lw=0.4, 
                              ax=ax[1], legend=None, colormap=cmap)
    sns.despine(fig=fig, top=False, left=True)
    
    for i in [0,1]:
        
        for x, bar in enumerate(ax[i].patches):
        
            if x in range(16,32):
                bar.set_hatch('////')
        
        ax[i].margins(x=0, y=0)
        handles, labels = ax[0].get_legend_handles_labels()
        ax[0].legend(reversed(handles), reversed(labels), markerscale=5, 
                     fontsize=16, bbox_to_anchor=[.4, -0.4, 0.85, 0.5])
        ax[i].set_xticklabels(tissue_type_long, rotation=45, fontsize=16)
        ax[i].set_ylabel('Percent', fontsize=16)
        ax[i].tick_params('y', labelsize=14)
        ax[0].set_title("Young", fontsize=18)
        ax[1].set_title("Old", fontsize=18)
        
    return fig


if __name__ == "__main__":
    data = pd.read_csv("data/Mouse_aging_mtDNA_summary.csv")
    
    young = process_data(data, 4.5)
    old = process_data(data, 26)
    
    fig = plot_figure()
    
    if not os.path.isdir('figures/'):
        os.mkdir('figures/')
        
    fig.savefig("figures/Supp_Figure_6.png", facecolor='white', dpi=600, 
                bbox_inches='tight') 
    fig.savefig("figures/Supp_Figure_6.pdf", facecolor='white', dpi=600, 
                bbox_inches='tight') 

    
    

