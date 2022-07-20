#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 09:18:01 2022

@author: scottrk
"""
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from GlobalVars_ import color_cycle, tissue_type, tissue_type_long, \
    mut_type, mut_type_pretty, mut_type_conv
from HelperFuncs_ import p_val_convert
from compile_data import melt_summary, summary_import
from compute_stats import heatmap_stats


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

    fig = plt.figure(layout='tight', figsize=(14, 18))
    top, bottom = fig.subfigures(nrows=2, ncols=1)
    axd = top.subplot_mosaic(mosaic, gridspec_kw={'hspace': 0.1})
    axd2 = bottom.subplot_mosaic(mosaic2, gridspec_kw={'hspace': 0.1})

    top.subplots_adjust(hspace=0.1)
    axd2ins = inset_axes(axd2['A'],
                         width="32.5%",
                         height="100%",
                         bbox_to_anchor=[-0.333, 1.5, 1, 1],
                         bbox_transform=axd2['A'].transAxes)

    return fig, axd, axd2, axd2ins


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

    fig, axd, axd2, axd2ins = setup_figure()

    spectrum1 = spectrum(data.query("Age=='Young'"), mut_type, axd['A'], 
                         ymin=0, ymax=5e-6)

    spectrum2 = spectrum(data.query("Age=='Old'"), mut_type, axd2['A'], 
                         fill=True, ymin=0, ymax=1.55e-5)

    sns.barplot(x='Class', y='Frequency', hue='Tissue',
                data=data.query("Age=='Old' & Class in ['C>A/G>T', 'C>G/G>C']"),
                order=mut_type[2:4], ci='sd', edgecolor='black', lw=1.5,
                errwidth=1.5, capsize=0.04, errcolor='black', ax=axd2ins)

    fig.canvas.draw()
    plt.setp(axd2ins.get_yaxis().get_offset_text(), visible=False)
    order_of_mag = axd2ins.get_yaxis().get_offset_text().get_text()[-2:]
    string = "Mutation Frequency ($\mathregular{10^{" + str(order_of_mag) + "}}$)"
    axd2ins.set_ylabel(string, fontsize=15)
    axd2ins.set_xticklabels(mut_type_pretty[2:4], fontsize=15)
    axd2ins.tick_params('y', labelsize=15)
    axd2ins.legend_.remove()
    axd2ins.set_xlabel('')

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

        heatmap = mod_heatmap(p_val_convert(young_tissue_comp), tissue_type[:-1],
                              axd[subplot[i]])

        heatmap2 = mod_heatmap(p_val_convert(old_tissue_comp), tissue_type[:-1],
                               axd2[subplot[i]])

        if subplot[i] != 'B':
            heatmap.set_yticklabels('')
            heatmap2.set_yticklabels('')
            heatmap.tick_params(labelrotation=0, labelsize=14)
            heatmap2.tick_params(labelrotation=0, labelsize=14)
        else:
            heatmap.tick_params(labelrotation=0, labelsize=14)
            heatmap2.tick_params(labelrotation=0, labelsize=14)

    for i, ax in enumerate(axd2ins.patches):
        j = i // 2
        axd2ins.patches[i].set_facecolor(color_cycle[j])
        axd2ins.patches[i].set_edgecolor('black')

    if not os.path.isdir("figures/"):
        os.mkdir("figures")

    fig.savefig('figures/Figure_2.png', dpi=600, facecolor='white')
