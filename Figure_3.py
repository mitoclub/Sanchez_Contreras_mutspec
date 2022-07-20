#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 11:50:07 2022

@author: scottrk
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as Patch
import seaborn as sns
from compute_stats import Fig3_stats
from compile_data import summary_import, melt_summary
from GlobalVars_ import tissue_type, tissue_type_long, mut_type, \
    mut_type_pretty, color_cycle, R_lib_path


def mod_ratios(mratios_file):
    ratio_df = pd.read_csv(mratios_file)
    ratio_df["Estimate"] = np.log2(ratio_df["Estimate"])
    ratio_df["confInt1"] = np.log2(ratio_df["confInt1"])
    ratio_df["confInt2"] = np.log2(ratio_df["confInt2"])
    ratio_df["Estimate"] = np.where(ratio_df['Class'] == "A>C/T>G", 0, ratio_df['Estimate'])
    ratio_df["confInt1"] = np.where(ratio_df['Class'] == "A>C/T>G", 0, ratio_df['confInt1'])
    ratio_df["confInt2"] = np.where(ratio_df['Class'] == "A>C/T>G", 0, ratio_df['confInt2'])
    ratio_df["qVal"] = np.where(ratio_df['Class'] == "A>C/T>G", np.nan, ratio_df['qVal'])
    ratio_df['qVal'] = np.where(ratio_df['qVal'] > 0.05, 0, ratio_df['qVal'])
    ratio_df['qVal'] = np.where((ratio_df['qVal'] > 0) & (ratio_df['qVal'] < 0.0001), 4, ratio_df['qVal'])
    ratio_df['qVal'] = np.where((0.0001 < ratio_df['qVal']) & (ratio_df['qVal'] < 0.001), 3, ratio_df['qVal'])
    ratio_df['qVal'] = np.where((0.001 < ratio_df['qVal']) & (ratio_df['qVal'] < 0.01), 2, ratio_df['qVal'])
    ratio_df['qVal'] = np.where((0.01 < ratio_df['qVal']) & (ratio_df['qVal'] < 0.05), 1, ratio_df['qVal'])
    ratio_df['u_confInt_delta'] = ratio_df['confInt2'] - ratio_df['Estimate']
    ratio_df['l_confInt_delta'] = ratio_df['Estimate'] - ratio_df['confInt1']

    return ratio_df


def setup_figure():
    mosaic = """AAAAAA
                AAAAAA
                AAAAAA
                AAAAAA
                BCDEFG
                .....G"""

    fig = plt.figure(figsize=(14, 8))
    axd = fig.subplot_mosaic(mosaic, gridspec_kw={'hspace': 1})

    axd['A'].set_ylim(-1.5, 3)
    axd['A'].axhline(0, color='black')
    axd['A'].margins(x=0.01)
    axd['A'].set_xticklabels(mut_type_pretty, fontsize=18)
    axd['A'].set_yticks([-2, -1, 0, 1, 2, 3])
    axd['A'].set_yticklabels([-4, -2, 1, 2, 4, 8], fontsize=18)

    return fig, axd


if __name__ == "__main__":

    if not os.path.isfile("data/stats/Figure_3_ratio_statistics.csv"):
        if not os.path.isfile("data/imported_data/summary_data_tidy.csv"):
            data = melt_summary(summary_import("data/Mouse_aging_mtDNA_summary.csv"))
            data.to_csv("data/imported_data/summary_data_tidy.csv")

        if not os.path.isdir("data/stats/"):
            os.mkdir("data/stats/")
        else:
            Fig3_stats(R_lib_path)
    else:
        ratio_df = mod_ratios("data/stats/Figure_3_ratio_statistics.csv")

    sort_dict = {'K': 0, 'L': 1, 'EC': 3, 'R': 4, 'Hi': 5, 'C': 6, 'M': 7, 'He': 8}
    ratio_df.sort_values(by=['Tissue'], key=lambda x: x.map(sort_dict),
                         inplace=True)

    fig, axd = setup_figure()

    sns.barplot(x='Class', y='Estimate', hue='Tissue', data=ratio_df,
                palette='bright', order=mut_type, ci='sd', edgecolor='black',
                lw=1.5, ax=axd['A'])

    sort_dict2 = {'C>T/G>A': 0, 'A>G/T>C': 1, 'C>A/G>T': 2, 'C>G/G>C': 3, 'A>T/T>A': 4, 'A>C/T>G': 5}
    resorted_df_list = []

    for tissue in tissue_type[:-1]:
        resorted_df = ratio_df.query("Tissue==@tissue").sort_values(by=['Class'],
                                                                    key=lambda x: x.map(sort_dict2))
        resorted_df_list.append(resorted_df)

    resorted_df = pd.concat(resorted_df_list)
    x_coords = []

    for i, ax in enumerate(axd['A'].patches):
        j = i // 6
        axd['A'].patches[i].set_facecolor(color_cycle[j])
        axd['A'].patches[i].set_edgecolor('black')
        legend = [Patch.Patch(facecolor=color_cycle[x], edgecolor='black',
                              label=tissue_type_long[x]) for x in range(8)]
        x_coords.append(ax.get_x() + 0.5 * ax.get_width())

    upper = list(resorted_df['u_confInt_delta'])
    lower = list(resorted_df['l_confInt_delta'])

    axd['A'].errorbar(x=x_coords, y=resorted_df['Estimate'], yerr=lower,
                      fmt="none", c="black", capsize=4)

    axd['A'].legend(handles=legend, ncol=4, fontsize=16, bbox_to_anchor=(.1, 1),
                    frameon=False)

    axd['A'].set_ylabel("Fold Change w/Age ($log_2$ Scaled)", fontsize=18)
    axd['A'].tick_params(labelsize=18)
    axd['A'].set_xlabel("")

    subplot = ['B', 'C', 'D', 'E', 'F']

    for i, mut_class in enumerate(mut_type[:-1]):
        sns.heatmap(pd.DataFrame(ratio_df.query('Class==@mut_class')['qVal']).T,
                    ax=axd[subplot[i]], square=True, vmin=0, vmax=10,
                    cbar=False, cmap='rocket_r', linewidths=0.5,
                    linecolor='black')

        sns.despine(top=False, right=False, ax=axd[subplot[i]])

        axd[subplot[i]].set_xticklabels(tissue_type[:-1], fontsize=14)
        axd[subplot[i]].set_yticklabels('')

    sns.heatmap(pd.DataFrame([4, 3, 2, 1, 0]), square=True, cbar=False,
                cmap='rocket_r', vmin=0, vmax=10,
                yticklabels=['q<0.0001', '0.0001<q<0.001', '0.001<q<0.01', '0.01<q<0.05', 'ns'],
                xticklabels=[], ax=axd["G"])

    axd['G'].tick_params(labelsize="large", labelleft=False, labelright=True,
                         left=False, right=True, labelrotation=0)
    axd['B'].set_yticklabels(["q-value"], rotation=0, fontsize=14)
    axd['A'].text(x=4.9, y=0.1, s="ND", fontsize=25)
    pos = axd['G'].get_position()
    pos.y0 = 0.08
    axd['G'].set_position(pos)

    if not os.path.isdir("figures"):
        os.mkdir("figures/")

    fig.savefig("figures/Figure_3.png", facecolor="white", dpi=600)
