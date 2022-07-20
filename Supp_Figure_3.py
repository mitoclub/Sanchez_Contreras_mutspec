#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 16:13:51 2022

@author: scottrk
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as Patch
from matplotlib.lines import Line2D
import seaborn as sns
from GlobalVars_ import young_mouse_id, old_mouse_id, tissue_type, \
    tissue_type_long, color_cycle
from compile_data import summary_import


def collate_copy_num_data(copy_num_table, summary_data):
    freq_SNV_data = []
    freq_InDel_data = []
    copy_num_data = []
    tissue_list = []

    for group in [young_mouse_id, old_mouse_id]:

        for mouse in group:
            mouseid = mouse.split('_')[0]

            for tissue in tissue_type[:-1]:

                if mouseid in list(copy_num_table['Mouse_ID']):
                    try:
                        copy_num_val = copy_num_table.query("Mouse_ID==@mouseid & Tissue_Abbrev==@tissue")[
                            'mt_tert_ratio']

                        if str(list(copy_num_val)[0]) != 'nan':
                            freq_SNV_data.append(
                                float(summary_data.query("MouseID==@mouseid & Tissue==@tissue")['Total_SNV_Freq']))
                            freq_InDel_data.append(
                                float(summary_data.query("MouseID==@mouseid & Tissue==@tissue")['Total_InDel_Freq']))
                            copy_num_data.append(float(copy_num_val))
                            tissue_list.append(tissue)
                    except:
                        pass

    df = pd.DataFrame([freq_SNV_data,
                       freq_InDel_data,
                       copy_num_data,
                       tissue_list],
                      index=["SNV_Freq", "InDel_Freq", "Copy_Num", "Tissue"]
                      ).T

    df.query("Tissue!='Blood'", inplace=True)

    return df


def setup_figure():
    fig, axs = plt.subplots(nrows=3, layout='constrained', figsize=(6, 12))

    axs[1].margins(x=0, y=0)
    axs[2].margins(x=0, y=0)

    return fig, axs


def plot_subplot_A(ax):
    sns.boxplot(x="Tissue_Abbrev", y="mt_tert_ratio", hue="Group",
                hue_order=['Young', 'Aged', 'ELAM'], data=copy_num_data,
                order=tissue_type[:-1], flierprops=dict(marker='o', color='black'),
                whiskerprops=dict(color='black'), medianprops=dict(color='black'),
                linewidth=1.5, ax=ax)

    ax.set_ylim(0, 15000)

    ax.set_xlabel("Tissue", fontsize=16)
    ax.set_ylabel("mtDNA:nDNA", fontsize=14)
    ax.set_xticklabels(tissue_type_long, fontdict={'rotation': 45}, fontsize=14)

    legend = [Patch.Patch(facecolor='white', edgecolor='black', label='Young'),
              Patch.Patch(facecolor='lightgrey', edgecolor='black', label='Old'),
              Patch.Patch(facecolor='lightgrey', edgecolor='black', label='Elam', hatch='///')]
    ax.legend(handles=legend, title='Group')

    fc = ['#FFFFFF', '#FFFFFF', '#FFFFFF', '#ff9200', '#FFFFFF', '#ff9200',
          '#FFFFFF', '#0433ff', '#0433ff', '#FFFFFF', '#ff84ff', '#ff84ff',
          '#FFFFFF', '#76d5ff', '#76d5ff', '#FFFFFF', '#932191', '#932191',
          '#FFFFFF', '#935200', '#935200', '#FFFFFF', '#008e00', '#008e00',
          '#FFFFFF', '#ff2600', '#ff2600'
          ]

    ec = ['black', '#ff9200', 'black', 'black', 'black', 'black', '#0433ff',
          'black', 'black', '#ff84ff', 'black', 'black', '#76d5ff', 'black',
          'black', '#932191', 'black', 'black', '#935200', 'black', 'black',
          '#008e00', 'black', 'black', '#ff2600', 'black', 'black'
          ]

    for i, bar in enumerate(ax.patches):

        ax.patches[i].set_facecolor(fc[i])
        ax.patches[i].set_edgecolor(ec[i])

        if i in [5, 8, 11, 14, 17, 20, 23, 26]:
            bar.set_hatch('///')


def plot_subplot_B_C(ax, y, ylabel, xlabel="mtDNA:nDNA"):
    sns.scatterplot(x='Copy_Num', y=y, hue='Tissue', data=collated_data,
                    hue_order=tissue_type[:-1], palette=color_cycle, ax=ax)

    ax.legend(handles=[Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=color_cycle[i],
                              markersize=10) for i in range(0, 8)],
              labels=tissue_type_long,
              title="Tissue",
              ncol=2)

    plt.setp(ax.get_yaxis().get_offset_text(), visible=False)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=14)


if __name__ == "__main__":
    copy_num_data = pd.read_csv("data/misc_items/mito_copy_number_data.csv")

    collated_data = collate_copy_num_data(copy_num_data,
                                          summary_import("data/Mouse_aging_mtDNA_summary.csv")
                                          )
    fig, axs = setup_figure()

    plot_subplot_A(axs[0])
    plot_subplot_B_C(axs[1], y="SNV_Freq",
                     ylabel="SNV Frequency ($\mathregular{10^{-6}}$)")

    plot_subplot_B_C(axs[2], y="SNV_Freq",
                     ylabel="In/Del Frequency ($\mathregular{10^{-6}}$)")

    if not os.path.isdir("figures"):
        os.mkdir("figures/")

    fig.savefig('figures/Supp_Figure_3.png', dpi=600, facecolor='white')
