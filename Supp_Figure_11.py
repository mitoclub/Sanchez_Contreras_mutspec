#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 11:26:59 2022

@author: scottrk
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as Patch
import seaborn as sns
from compile_data import summary_import, melt_summary
from GlobalVars_ import color_cycle, tissue_type, tissue_type_long

if not os.path.isfile("data/imported_data/summary_data_tidy.csv"):
    if not os.path.isdir("data/imported_data/"):
        os.mkdir("data/imported_data")

    data = melt_summary(summary_import("data/Mouse_aging_mtDNA_summary.csv"))
    data.to_csv("data/imported_data/summary_data_tidy.csv")
else:
    data = pd.read_csv("data/imported_data/summary_data_tidy.csv", index_col=0)

data = data.query("Treatment!='perf' & Age=='Old' & Class=='Total_SNV_Freq'")

fig, ax = plt.subplots(ncols=1, figsize=(12, 7))

sns.barplot(x="Tissue", y="Frequency", hue="Treatment", data=data,
            order=tissue_type[:-1], palette='bright', ci='sd', edgecolor='black',
            lw=1.2, errwidth=1.5, capsize=0.1, errcolor='black', ax=ax)

hatches = ['///', '---']
fig_color = color_cycle * 3

for x, bar in enumerate(ax.patches):

    ax.patches[x].set_facecolor(fig_color[x])

    if x in range(8, 16):
        bar.set_hatch(hatches[0])
    elif x in range(16, 23):
        bar.set_hatch(hatches[1])

patch1 = Patch.Patch(facecolor='white', edgecolor='black', label='Old')

patch2 = Patch.Patch(facecolor='white', edgecolor='black', hatch=hatches[0],
                     label='Old+Elam.')

patch3 = Patch.Patch(facecolor='white', edgecolor='black', hatch=hatches[1],
                     label='Old+NMN')

legend = [patch1, patch2, patch3]

ax.legend(handles=legend, markerscale=5, fontsize=16)

ax.set_xticklabels(tissue_type_long, fontsize=14)
ax.set_xlabel("Tissue", fontsize=18)
plt.setp(ax.get_yaxis().get_offset_text(), visible=False)
ax.set_ylabel("Mutation Frequency ($\mathregular{10^{-6}}$)", fontsize=18)
ax.tick_params('y', labelsize=14)

if not os.path.isdir("figures"):
    os.mkdir("figures/")

fig.savefig('figures/Supp_Figure_11.png', dpi=600, facecolor='white')
