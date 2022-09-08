#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 15:19:11 2022

@author: scottrk
"""
import os
import pandas as pd
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as Patch
import seaborn as sns
from GlobalVars_ import old_mouse_id, perfused_mouse_id
from compile_data import get_depth_data


def get_clone_count(df, alt_count):
    clone_count = df.query(" alt_count>@alt_count & ref!='-' & alt!='-' ")
    denovo_count = df.query(" alt_count<@alt_count & ref!='-' & alt!='-' ")
    total_vars = df.query(" ref!='-' & alt!='-' ")

    return len(clone_count), len(denovo_count), len(total_vars)


cohort_dict = {0: 'Old', 1: 'Perf'}
data_list = []

for index, group in enumerate([old_mouse_id, perfused_mouse_id]):

    for tissue in ['K', 'L', 'Hi', 'C', 'M', 'B']:

        for mouse in group:

            try:
                mut_data = pd.read_csv(glob.glob('data/mut_files/' + mouse + '_' + tissue + '*.mut')[0], sep='\t')
                mean_depth = get_depth_data(glob.glob('data/depth_files/' + mouse + '_' + tissue + '_' + '*.dcs'
                                                                                                         '.region'
                                                                                                         '.mutpos'
                                                                                                         '.vcf_depth'
                                                                                                         '.txt')[0])
                clone_count, denovo_count, total_var_count = get_clone_count(mut_data, 2)
                data_list.append(pd.DataFrame([mouse, tissue, cohort_dict[index],
                                               clone_count, total_var_count,
                                               mean_depth, clone_count / mean_depth,
                                               clone_count / total_var_count]).T)
            except:
                pass

clone_freq_df = pd.concat(data_list, ignore_index=True)
clone_freq_df.columns = ["MouseID", "Tissue", "Cohort", "Clone_Count",
                         "Total_Var_Count", "Mean_Depth", "Clone_Freq",
                         "Clone_Percentage"]

fig, ax = plt.subplots(ncols=1, constrained_layout=True, figsize=(15, 8))

color_cycle = ['#ff9200', '#0433ff', '#932191', '#935200', '#008e00', 'gold']
fc = color_cycle * 2

plot = sns.barplot(x='Tissue', y='Clone_Freq', hue='Cohort', data=clone_freq_df,
            ci='sd', lw=1.2, ec='black', errwidth=1.5, capsize=0.1, errcolor='black')

sns.stripplot(x="Tissue", y="Clone_Freq", hue="Cohort", data=clone_freq_df,
              order=['K', 'L', 'Hi', 'C', 'M', 'B'], dodge=True, ax=ax,
              alpha=0.7, color='black')

for i in range(12):

    plot.patches[i].set_facecolor(fc[i])
    #plot.patches[i].set_edgecolor(ec[i])

    if i < 5:
        r, g, b, a = plot.patches[i].get_facecolor()
        plot.patches[i].set_facecolor((r, g, b, .15))  
        plot.patches[i].set_edgecolor((r, g, b, 1))
        plot.patches[i].set(lw=2.4)

    legend = [Patch.Patch(facecolor='lightgrey', edgecolor='black', label='Old'),
              Patch.Patch(facecolor='dimgrey', edgecolor='black', label='Old Perfused')]
    plot.legend(handles=legend, fontsize='xx-large', loc='upper right')
    
ax.set_ylabel("Clone Frequency", fontsize=22)
ax.set_xlabel("Tissue", fontsize=22)
ax.set_xticklabels(['Kidney', 'Liver', 'Hippo.', 'Cerebel.', 'Sk. Muscle', 'Blood'],
                   fontsize=18)
ax.tick_params('y', labelsize=18)

if not os.path.isdir("figures"):
    os.mkdir("figures/")

fig.savefig("figures/Supp_Figure_7.png", facecolor='white', dpi=600)
fig.savefig("figures/Supp_Figure_7.pdf", facecolor='white', dpi=600)
