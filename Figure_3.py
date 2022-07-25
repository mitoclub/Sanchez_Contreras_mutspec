#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 08:28:07 2022

@author: scottrk
"""
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from dna_features_viewer import BiopythonTranslator
import numpy as np
import pandas as pd
from operator import itemgetter
import os
from GlobalVars_ import tissue_type, tissue_type_long, color_cycle
from compile_data import mut_file_import, calc_clone_numbers


class MyCustomTranslator(BiopythonTranslator):

    def compute_feature_color(self, feature):

        if feature.type == "gene":
            return "purple"
        elif feature.type == "tRNA":
            return "blue"
        elif feature.type == "rRNA":
            return "orange"
        else:
            return "green"

    def compute_feature_label(self, feature):

        if feature.type is not None:
            return None
        # elif feature.type == "CDS":
        #    return "CDS here"
        else:
            return BiopythonTranslator.compute_feature_label(self, feature)

    def compute_filtered_features(self, features):
        """Do not display promoters. Just because."""
        return [
            feature for feature in features
            if feature.type is not None
        ]


def setup_figure():
    mosaic = """AABB
                CCDD
                EEFF
                GGHH
                IIJJ"""

    fig = plt.figure(figsize=(12, 14), constrained_layout=True, facecolor='white')
    axd = fig.subplot_mosaic(mosaic)

    fig.supylabel('Variant Allele Fraction (Log$_2$ Scaling)', y=0.35, size=20)
    fig.supxlabel('Genome Position', size=20, x=0.55)
    axdins = inset_axes(axd['A'],
                        width="32.5%",
                        height="100%",
                        bbox_to_anchor=[-0.095, .5, 1.1, .5],
                        bbox_transform=axd['A'].transAxes)

    return fig, axd, axdins


def plot_clone_numbers(x, y, hue, data, order, ax, ylim, xticklabels, xlabel, ylabel, fc, ec, legend=True):
    sns.barplot(x=x, y=y, hue=hue, data=data, order=order, ci='sd',
                       lw=1.2, errwidth=1.5, capsize=0.1, errcolor='black',
                       ax=ax)
    
    sns.despine(ax=ax)
    
    ax.set_xticklabels(xticklabels, rotation=45, fontsize=14)
    ax.tick_params('y', labelsize=14)
    ax.set_xlabel(xlabel, fontsize=16)
    plt.setp(ax.get_yaxis().get_offset_text(), visible=False)
    ax.set_ylabel(ylabel, fontsize=16)

    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])

    for i in range(len(fc)):
        ax.patches[i].set_facecolor(fc[i])
        ax.patches[i].set_edgecolor(ec[i])

        if i < len(fc) / 2:
            ax.patches[i].set(lw=2.4)

    if legend:
        new_legend = [Patch.Patch(facecolor='white', edgecolor='black', label='Young'),
                      Patch.Patch(facecolor='lightgrey', edgecolor='black', label='Old')]
        ax.legend(handles=new_legend, fontsize='x-large', frameon=False, ncol=2, bbox_to_anchor=[0.5, 0.73, 0.5, 0.5])
    elif not legend:
        ax.legend_.remove()
    else:
        pass


def clone_plot(df, tissue, flat_file, ax):
    graphic_record = MyCustomTranslator().translate_record(flat_file)
    graphic_record.feature_level_height = 0
    graphic_record.plot(figure_width=1, ax=ax)
    y = ax.twinx()

    for index, group in enumerate(['Young', 'Old']):
        data = df.query("Cohort==@group & Tissue==@tissue & \
                        alt_count>2 & alt!='-' & ref!='-' & VAF<0.01")

        data['VAF'] = np.log2(data['VAF'].mul(10000).astype(float))

        if index == 0:
            bottom = 0.67

        if index == 1:
            data['VAF'] = data['VAF'].mul(-1)
            bottom = -0.67

        markerline, stemlines, baseline = y.stem(data['start'], data['VAF'],
                                                 linefmt='C' + str(index) + '-',
                                                 markerfmt='C' + str(index) + 'o',
                                                 basefmt=" ", label=group,
                                                 bottom=bottom)
        stemlines.set_linewidths(1)
        markerline.set_markersize(3.4)

    for x in [1, 2, 3, 4, 5, 6]:
        y.axhline(y=-x, c='black', linestyle='--', alpha=0.11, zorder=0)
        y.axhline(y=x, c='black', linestyle='--', alpha=0.11, zorder=0)

    y.axhline(y=0, c='black')

    y.tick_params(labelsize='x-large')
    y.spines['bottom'].set_color('black')
    y.set_yticks([-6, 0, 6])

    if i % 2 == 1:
        y.set_yticklabels([])

    else:
        y.set_yticklabels(['0.0064', '0', '0.0064'])

    ax.set_yticks([1, 0, 1])
    ax.set_ylim(-1, 1)
    ax.set_xticklabels(['0', '2k', '4k', '6k', '8k', '10k', '12k', '14k', '16k'],
                       fontsize=16, color='black')
    y.spines["left"].set_position(("axes", 0))
    y.yaxis.set_label_position('left')
    y.yaxis.set_ticks_position('left')

    return ax


if __name__ == "__main__":

    if not os.path.isfile("data/imported_data/summary_clone_data.csv"):
        if not os.path.isdir("data/imported_data"):
            os.mkdir("data/imported_data/")

        mut_data = mut_file_import()
        final_clone_data = calc_clone_numbers(mut_data)

        mut_data.to_csv("data/imported_data/mut_file_data.csv")
        final_clone_data.to_csv("data/imported_data/summary_clone_data.csv")
    else:
        mut_data = pd.read_csv("data/imported_data/mut_file_data.csv",
                               index_col=[0])
        final_clone_data = pd.read_csv("data/imported_data/summary_clone_data.csv")

    final_clone_data = final_clone_data.query("Cohort in ['Young', 'Old'] & Tissue != 'B'")
    

    fig, axd, axdins = setup_figure()

    fc = ['white'] * 8 + color_cycle
    ec = color_cycle + ['black'] * 8

    plot_clone_numbers(x="Tissue", y="Clone_Freq", hue="Cohort",
                       data=final_clone_data, order=tissue_type[:-1], ylim=None,
                       ax=axd['A'], xlabel='',
                       ylabel="Clone Frequency($\mathregular{10^{-2}}$)",
                       xticklabels=tissue_type_long, fc=fc, ec=ec)

    plot_clone_numbers(x="Tissue", y="Percent_Clone", hue="Cohort",
                       data=final_clone_data, order=tissue_type[:-1], ylim=[0, 14],
                       ax=axd['B'], xlabel='', ylabel="Percent Clones",
                       xticklabels=tissue_type_long, fc=fc, ec=ec)

    sub_data = final_clone_data.query("Tissue in ['C', 'M', 'He']")

    plot_clone_numbers(x="Tissue", y="Clone_Freq", hue="Cohort", data=sub_data,
                       order=tissue_type[-4:-1], ylim=None, ax=axdins,
                       xlabel='', ylabel="", xticklabels=['', '', ''],
                       fc=itemgetter(0, 1, 2, 13, 14, 15)(fc), ec=ec[5:11],
                       legend=False)
    
    sns.despine(ax=axdins, top=False, right=False)
    
    sns.set_palette("tab10")
    
    subplot = ['C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']

    for i, tissue in enumerate(tissue_type[:-1]):
        clone_plot(mut_data, tissue, "data/misc_items/NC_005089.1[1..16299].flat",
                   axd[subplot[i]])

        axd[subplot[i]].set_title(tissue_type_long[i], fontsize=20, y=1,
                                  backgroundcolor='white')

    sns.set_palette('bright')

    if not os.path.isdir("figures"):
        os.mkdir("figures/")

    fig.savefig("figures/Figure_3.png", facecolor="white", dpi=600)
