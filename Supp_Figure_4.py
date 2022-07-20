#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 17:17:18 2022

@author: scottrk
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as Patch
import seaborn as sns
from compile_data import summary_import

fc = ['#ffffff', '#ffffff', '#ffffff', '#ffffff', '#ffffff', '#ffffff',
      '#ff9200', '#0433ff', '#932191', '#935200', '#008e00', 'gold']

ec = ['#ff9200', '#0433ff', '#932191', '#935200', '#008e00', 'gold', '#000000',
      '#000000', '#000000', '#000000', '#000000', '#000000']

fig, ax = plt.subplots(ncols=1, constrained_layout=True, figsize=(15, 8))

if not os.path.isfile("data/imported_data/summary_data_wide.csv"):
    if not os.path.isdir("data/imported_data/"):
        os.mkdir("data/imported_data")

    summary_data = summary_import("data/Mouse_aging_mtDNA_summary.csv")
else:
    summary_data = pd.read_csv("data/imported_data/summary_data_wide.csv",
                               index_col=0)

data = summary_data.query("Treatment in ['NT', 'perf'] & Age=='Old' & Tissue not in ['EC', 'R', 'He']")

plot = sns.barplot(x="Tissue", y="Total_SNV_Freq", hue="Treatment", data=data,
                   order=['K', 'L', 'Hi', 'C', 'M', 'B'], palette='bright',
                   ci='sd', edgecolor='black', lw=1.2, errwidth=1.5,
                   capsize=0.1, errcolor='black', ax=ax)

sns.stripplot(x="Tissue", y="Total_SNV_Freq", hue="Treatment", data=data,
              order=['K', 'L', 'Hi', 'C', 'M', 'B'], dodge=True, ax=plot,
              alpha=0.7, color='black')

for i in range(12):

    plot.patches[i].set_facecolor(fc[i])
    plot.patches[i].set_edgecolor(ec[i])

    if i < 5:
        plot.patches[i].set(lw=2.4)

    legend = [Patch.Patch(facecolor='white', edgecolor='black', label='Old'),
              Patch.Patch(facecolor='lightgrey', edgecolor='black', label='Old Perfused')]
    plot.legend(handles=legend, fontsize='xx-large', loc='upper right')

plot.set_ylabel("SNV Mutation Frequency ($\mathregular{10^{-6}}$)", fontsize=20)
plot.set_xlabel("Tissue", fontsize=20)
plot.set_xticklabels(['Kidney', 'Liver', 'Hippo.', 'Cerebel.', 'Sk. Muscle', 'Blood'])
plot.tick_params(labelsize=18)

if not os.path.isdir("figures"):
    os.mkdir("figures/")

fig.savefig('figures/Supp_Figure_4.png', dpi=600, facecolor='white')
