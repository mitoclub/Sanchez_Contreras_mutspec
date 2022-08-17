#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 09:18:01 2022

@author: scottrk
"""
import os
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from GlobalVars_ import tissue_type, tissue_type_long, mut_type, \
     mut_type_pretty, color_cycle, mut_type_conv, R_lib_path, tissue_type_abbrev
from HelperFuncs_ import p_val_convert
from compile_data import melt_summary, summary_import
from compute_stats import heatmap_stats, Fig2C_stats


def spectrum(data, mut_type, ax, fill=False, ymin=None, ymax=None, legend=True):
    
    plot = sns.barplot(x='Class', y='Frequency', hue='Tissue', data=data,
                       palette='bright', order=mut_type, ci='sd',
                       edgecolor='black', lw=1.5, errwidth=1.5,
                       capsize=0.04, errcolor='black', ax=ax)

    sns.stripplot(x="Class", y="Frequency", hue="Tissue", data=data,
                  order=mut_type, ax=plot, alpha=0.7, dodge=True, color='black')

    sns.despine(ax=ax)

    if ymin is not None and ymax is not None:
        plot.set_ylim(ymin, ymax)

    if not fill:
        for i, ax in enumerate(plot.patches):

            j = i // 6
            plot.patches[i].set_facecolor('white')
            plot.patches[i].set_edgecolor(color_cycle[j])
            plot.patches[i].set(lw=2.4)

            if legend:
                legend = [Patch.Patch(facecolor='white', edgecolor=color_cycle[x],
                                      label=tissue_type_long[x]) for x in range(8)]
    else:

        for i, ax in enumerate(plot.patches):

            j = i // 6
            plot.patches[i].set_facecolor(color_cycle[j])
            plot.patches[i].set_edgecolor('black')

            if legend:
                legend = [Patch.Patch(facecolor=color_cycle[x], edgecolor='black',
                                      label=tissue_type_long[x]) for x in range(8)]

    fig.canvas.draw()
    plt.setp(plot.get_yaxis().get_offset_text(), visible=False)
    order_of_mag = plot.get_yaxis().get_offset_text().get_text()[-2:]
    string = "Mutation Frequency ($\mathregular{10^{" + str(order_of_mag) + "}}$)"
    plot.set_xticklabels(mut_type_pretty)
    plot.set_xlabel('')
    plot.set_ylabel(string, fontsize=20)
    plot.tick_params(labelsize=18)
    plot.margins(x=.01)
    plot.legend(handles=legend, ncol=4, fontsize=16, bbox_to_anchor=(.9, 1.2),
                frameon=False)

    return plot


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


def mod_heatmap(data, labels, ax):
    plot = sns.heatmap(data, square=True, cbar=False, cmap='rocket_r',
                       linewidths=.5, linecolor='black', vmin=0,
                       vmax=10, xticklabels=labels, yticklabels=labels, ax=ax)

    sns.despine(left=True, bottom=True, ax=plot)
    plot.set_facecolor('silver')
    plot.tick_params(labelsize=18)

    return plot


def setup_figure():
    mosaic = """AAAAAA
                AAAAAA
                BCDEFG"""

    mosaic2 = """AAAAAA
                 AAAAAA
                 BCDEFG"""
                 
    mosaic3 = """AAAAAA
                 AAAAAA
                 AAAAAA
                 AAAAAA
                 BCDEFG
                 .....G"""

    fig = plt.figure(layout='tight', figsize=(14, 22))
    top, middle, bottom = fig.subfigures(nrows=3, ncols=1)
    axd = top.subplot_mosaic(mosaic, gridspec_kw={'wspace':0.016, 'hspace': 0.3})
    axd2 = middle.subplot_mosaic(mosaic2, gridspec_kw={'wspace':0.016, 'hspace': 0.3})
    axd3 = bottom.subplot_mosaic(mosaic3, gridspec_kw={'wspace':0.08, 'hspace': 0.3})
    top.subplots_adjust(hspace=0.1)
    axd2ins = inset_axes(axd2['A'], width="32.5%", height="130%",
                         bbox_to_anchor=[-0.333, 1, 1, 1],
                         bbox_transform=axd2['A'].transAxes)
    axd3['A'].set_ylim(-1.5, 3)
    axd3['A'].axhline(0, color='black')
    axd3['A'].margins(x=0.01)
    axd3['A'].set_xticklabels(mut_type_pretty, fontsize=18)
    axd3['A'].set_yticks([-2, -1, 0, 1, 2, 3])
    axd3['A'].set_yticklabels([-4, -2, 1, 2, 4, 8], fontsize=18)
    
    return fig, axd, axd2, axd3, axd2ins

def Fig_2C(data, ax):
    plot = sns.barplot(x='Class', y='Estimate', hue='Tissue', data=ratio_df,
                palette='bright', order=mut_type, ci='sd', edgecolor='black',
                lw=1.5, ax=ax)
    
    sort_dict2 = {'C>T/G>A': 0, 'A>G/T>C': 1, 'C>A/G>T': 2, 'C>G/G>C': 3, 'A>T/T>A': 4, 'A>C/T>G': 5}
    resorted_df_list = []
    
    for tissue in tissue_type[:-1]:
        resorted_df = data.query("Tissue==@tissue").sort_values(by=['Class'],
                                                                    key=lambda x: x.map(sort_dict2))
        resorted_df_list.append(resorted_df)
    
    resorted_df = pd.concat(resorted_df_list)
    x_coords = []
    
    for i, patch in enumerate(ax.patches):
        j = i // 6
        patch.set_facecolor(color_cycle[j])
        patch.set_edgecolor('black')
        x_coords.append(patch.get_x() + 0.5 * patch.get_width())
    legend = [Patch.Patch(facecolor=color_cycle[x], edgecolor='black',
                          label=tissue_type_long[x]) for x in range(8)]

    
    upper = list(resorted_df['u_confInt_delta'])
    lower = list(resorted_df['l_confInt_delta'])
    
    plot.errorbar(x=x_coords, y=resorted_df['Estimate'], yerr=lower,
                      fmt="none", c="black", capsize=4)
    
    plot.legend(handles=legend, ncol=4, fontsize=16, bbox_to_anchor=(.1, 1),
                    frameon=False)
    
    plot.set_ylabel("Fold Change w/Age ($log_2$ Scaled)", fontsize=18)
    plot.tick_params(labelsize=18)
    plot.set_xlabel("")

    return plot

def Fig_2C_heatmap(data, mut_class, ax):

    plot = sns.heatmap(pd.DataFrame(data.query('Class==@mut_class')['qVal']).T,
                ax=ax, square=True, vmin=0, vmax=10,
                cbar=False, cmap='rocket_r', linewidths=0.5,
                linecolor='black')

    sns.despine(top=False, right=False, ax=axd[subplot[i]])

    plot.set_xticklabels(tissue_type_abbrev[:-1], fontsize=14)
    plot.set_yticklabels('')
        
    return plot

if __name__ == '__main__':

    custom_dict = {'K': 0, 'L': 1, 'EC': 3, 'R': 4, 'Hi': 5, 'C': 6, 'M': 7, 'He': 8}

    if not os.path.isfile("data/imported_data/summary_data_tidy.csv"):
        if not os.path.isdir("data/imported_data/"):
            os.mkdir("data/imported_data")

        data = melt_summary(summary_import("data/Mouse_aging_mtDNA_summary.csv"))
        data.to_csv("data/imported_data/summary_data_tidy.csv")
    else:
        data = pd.read_csv("data/imported_data/summary_data_tidy.csv", index_col=0)
        data = data.query("Treatment=='NT'").sort_values(by=['Tissue'],
                                                         key=lambda x: x.map(custom_dict))

    if not os.path.isfile("data/stats/Figure_2_ratio_statistics.csv"):

        if not os.path.isdir("data/stats/"):
            os.mkdir("data/stats/")
        else:
            Fig2C_stats(R_lib_path)
            ratio_df = mod_ratios("data/stats/Figure_2_ratio_statistics.csv")
    else:
        ratio_df = mod_ratios("data/stats/Figure_2_ratio_statistics.csv")

    sort_dict = {'K': 0, 'L': 1, 'EC': 3, 'R': 4, 'Hi': 5, 'C': 6, 'M': 7, 'He': 8}
    ratio_df.sort_values(by=['Tissue'], key=lambda x: x.map(sort_dict),
                         inplace=True)

    fig, axd, axd2, axd3, axd2ins = setup_figure()

    spectrum1 = spectrum(data.query("Age=='Young'"), mut_type, axd['A'], 
                         ymin=0, ymax=5e-6)

    spectrum2 = spectrum(data.query("Age=='Old'"), mut_type, axd2['A'], 
                         fill=True, ymin=0, ymax=1.55e-5)
    
    sns.barplot(x='Class', y='Frequency', hue='Tissue',
                data=data.query("Age=='Old' & Class in ['C>A/G>T', 'C>G/G>C']"),
                order=mut_type[2:4], ci='sd', edgecolor='black', lw=1.5,
                errwidth=1.5, capsize=0.04, errcolor='black', ax=axd2ins)

    spectrum_change = Fig_2C(ratio_df, axd3['A'])
    
    fig.canvas.draw()
    plt.setp(axd2ins.get_yaxis().get_offset_text(), visible=False)
    order_of_mag = axd2ins.get_yaxis().get_offset_text().get_text()[-2:]
    string = "Mutation Freq. ($\mathregular{10^{" + str(order_of_mag) + "}}$)"
    axd2ins.set_ylabel(string, fontsize=15)
    axd2ins.set_xticklabels(mut_type_pretty[2:4], fontsize=15)
    axd2ins.tick_params('y', labelsize=15)
    axd2ins.legend_.remove()
    axd2ins.set_xlabel('')
    
    for i, ax in enumerate(axd2ins.patches):
        j = i // 2
        axd2ins.patches[i].set_facecolor(color_cycle[j])
        axd2ins.patches[i].set_edgecolor('black')
    

    subplot = ['B', 'C', 'D', 'E', 'F', 'G']

    for i, mut_class in enumerate(mut_type):

        age = 'Young'

        if not os.path.isfile("data/stats/Figure_2_" + age + "_" + mut_type_conv[mut_class] + "_stats.csv"):
            if not os.path.isdir("data/stats/"):
                os.mkdir("data/stats/")

            young_tissue_comp = heatmap_stats(data, age, mut_class,
                                              "data/stats/Figure_2_" + age + "_" + mut_type_conv[
                                                  mut_class] + "_stats.csv")
        else:
            young_tissue_comp = pd.read_csv(
                "data/stats/Figure_2_" + age + "_" + mut_type_conv[mut_class] + "_stats.csv",
                index_col=0)

        age = 'Old'

        if not os.path.isfile("data/stats/Figure_2_" + age + "_" + mut_type_conv[mut_class] + "_stats.csv"):
            if not os.path.isdir("data/stats/"):
                os.mkdir("data/stats/", mode=0o666)

            old_tissue_comp = heatmap_stats(data, age, mut_class,
                                            "data/stats/Figure_2_" + age + "_" + mut_type_conv[
                                                mut_class] + "_stats.csv")
        else:
            old_tissue_comp = pd.read_csv("data/stats/Figure_2_" + age + "_" + mut_type_conv[mut_class] + "_stats.csv",
                                          index_col=0)

        heatmap = mod_heatmap(p_val_convert(young_tissue_comp), tissue_type_abbrev[:-1],
                              axd[subplot[i]])

        heatmap2 = mod_heatmap(p_val_convert(old_tissue_comp), tissue_type_abbrev[:-1],
                               axd2[subplot[i]])

        heatmap3 = Fig_2C_heatmap(ratio_df, mut_class, axd3[subplot[i]])
        
        if subplot[i] != 'B':
            heatmap.set_yticklabels('')
            heatmap2.set_yticklabels('')
            heatmap.tick_params(labelrotation=0, labelsize=14)
            heatmap2.tick_params(labelrotation=0, labelsize=14)
        else:
            heatmap.tick_params(labelrotation=0, labelsize=14)
            heatmap2.tick_params(labelrotation=0, labelsize=14)


    sns.heatmap(pd.DataFrame([4, 3, 2, 1, 0]), square=True, cbar=False,
                cmap='rocket_r', vmin=0, vmax=10,
                yticklabels=['q<0.0001', '0.0001<q<0.001', '0.001<q<0.01', '0.01<q<0.05', 'ns'],
                xticklabels=[], ax=axd3["G"])
    axd3['G'].tick_params(labelsize="large", labelleft=False, labelright=True,
                         left=False, right=True, labelrotation=0)
    axd3['B'].set_yticklabels(["q-value"], rotation=0, fontsize=14)
    axd3['A'].text(x=4.9, y=0.1, s="ND", fontsize=25)
    pos = axd3['G'].get_position()
    pos.y0 = 0.08
    axd3['G'].set_position(pos)

    if not os.path.isdir("figures/"):
        os.mkdir("figures")

    fig.savefig('figures/Figure_2.png', dpi=600, facecolor='white')
    fig.savefig('figures/Figure_2.pdf', dpi=600, facecolor='white')
