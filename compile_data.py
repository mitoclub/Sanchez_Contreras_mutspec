#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 13:55:25 2022

@author: scottrk
"""
import pandas as pd
import glob
import os
import itertools
from collections import defaultdict
from GlobalVars_ import young_mouse_id, old_mouse_id, old_ss31_id, \
                        perfused_mouse_id, old_nmn_id, tissue_type, cohort


def summary_import(summary_file):
    """function takes the summary data file as output by the Duplex-Seq pipeline
    and converts it to a Pandas DataFrame for use in downstream figure generation."""

    summary_data = pd.read_csv(summary_file, usecols=['MouseID', 'Age', 'Treatment',
                                                      'Tissue', 'Total_SNV_Freq',
                                                      'Total_InDel_Freq',
                                                      'A>T/T>A_Freq', 'A>C/T>G_Freq',
                                                      'A>G/T>C_Freq', 'C>A/G>T_Freq',
                                                      'C>T/G>A_Freq', 'C>G/G>C_Freq'])

    summary_data.Age.replace([4.5, 26], ["Young", "Old"], inplace=True)

    summary_data.columns = ['MouseID', 'Age', 'Treatment', 'Tissue',
                            'Total_SNV_Freq', 'Total_InDel_Freq', 'A>T/T>A',
                            'A>C/T>G', 'A>G/T>C', 'C>A/G>T', 'C>T/G>A', 'C>G/G>C']

    return summary_data


def melt_summary(summary_data):
    summary_data_long = pd.melt(summary_data,
                                id_vars=['MouseID', 'Tissue', 'Treatment',
                                         'Age'], value_vars=['Total_SNV_Freq',
                                                             'Total_InDel_Freq', 'A>T/T>A',
                                                             'A>C/T>G', 'A>G/T>C', 'C>A/G>T',
                                                             'C>T/G>A', 'C>G/G>C'])

    summary_data_long.columns = ["MouseID", "Tissue", 'Treatment', 'Age', 'Class',
                                 'Frequency']

    return summary_data_long


def mut_file_import():
#    tissue_dict = {}
#    grouping_dict = {}
    data_list = []
    for index, group in enumerate([young_mouse_id, old_mouse_id, old_ss31_id, 
                                   old_nmn_id, perfused_mouse_id]):

        for tissue in tissue_type:
            
            for mouse in group:

                try:
                    
                    data = pd.read_csv(glob.glob('data/mut_files/' + mouse + '_' 
                                                 + tissue + '*.mut')[0], sep='\t')
                    
                    df = pd.concat([data, pd.Series([tissue]*len(data), name='Tissue')], 
                                   axis=1)
                    df = pd.concat([df, pd.Series([cohort[index]]*len(data), name='Cohort')], 
                                   axis=1)
                    data_list.append(df)

                except:
                    pass

    final_mut_df = pd.concat(data_list, axis=0)

    return final_mut_df


def get_depth_data(depth_file):
    df = pd.read_csv(depth_file, sep='\t', usecols=[1, 3])

    return df['DP'].median()


def calc_depth_plot_data(mouse_list, tissue, output=None):
    df_list = []

    temp_df = pd.DataFrame(index=[x for x in range(1, 16299)])

    for mouse in mouse_list:

        try:
            df = pd.read_csv(glob.glob('data/depth_files/'
                                       + mouse + '_' + tissue +
                                       '*.dcs.region.mutpos.vcf_depth.txt')[0],
                             sep='\t',
                             usecols=[1, 3])
            df.set_index("Pos", inplace=True)
            df.columns = [mouse]
            df_list.append(pd.concat([temp_df, df], axis=1))

        except:
            pass
    df = pd.concat(df_list, axis=1)
    df['Mean'] = df.mean(axis=1).round()
    df['Std_dev_upper'] = df['Mean'] + df.std(axis=1)
    df['Std_dev_lower'] = df['Mean'] - df.std(axis=1)

    if output is not None:
        df.to_csv(output)

    return df


def get_clone_count(df, alt_count, max_vaf):
    clone_count = len(df.query(" alt_count>@alt_count & ref!='-' & alt!='-' & VAF<=@max_vaf"))
    denovo_count = len(df.query(" alt_count<@alt_count & ref!='-' & alt!='-' "))
    total_vars = len(df.query("ref!='-' & alt!='-' & VAF<=@max_vaf"))

    clone_type_count = defaultdict(lambda: 0)

    for elmt in itertools.permutations(['A', 'T', 'C', 'G']):
        clone_type_count[elmt[0] + elmt[1]] = len(
            df.query("alt_count>@alt_count & ref==@elmt[0] & alt==@elmt[1] & VAF<=@max_vaf"))

        combined_clone_type = {'AT_TA': clone_type_count['AT'] + clone_type_count['TA'],
                               'AC_TG': clone_type_count['AC'] + clone_type_count['TG'],
                               'AG_TC': clone_type_count['AG'] + clone_type_count['TC'],
                               'GA_CT': clone_type_count['GA'] + clone_type_count['CT'],
                               'GT_CA': clone_type_count['GT'] + clone_type_count['CA'],
                               'GC_CG': clone_type_count['GC'] + clone_type_count['CG']
                               }

        combined_clone_type = pd.DataFrame(combined_clone_type, index=[0])

    return clone_count, denovo_count, total_vars, combined_clone_type


def calc_clone_numbers(mut_file_data):
    final_clone_data = pd.DataFrame(columns=["Mouse_ID", "Cohort", "Tissue", 
                                             "Clone_Count", "A>T/T>A_Count", 
                                             "A>C/T>G_Count", "A>G/T>C_Count",
                                             "G>A/C>T_Count", "G>T/C>A_Count", 
                                             "G>C/C>G_Count", "Percent_Clone", 
                                             "Clone_Freq"])

    for i, group in enumerate(cohort[:2]+cohort[-1:]):

        for mouse in [young_mouse_id, old_mouse_id, perfused_mouse_id][i]:
            
            
             
            for tissue in tissue_type:
                
                tissue_string = '_' + tissue + '_'
                
                try:

                    df = mut_file_data.query("sample.str.contains(@mouse.upper()) & \
                                             sample.str.contains(@tissue_string.upper())")
                    clone_count, denovo_count, total_vars, clone_type_count = get_clone_count(df, 2, 0.01)
                    percent_clone = 100 * (clone_count / total_vars)
    
                    median_depth = get_depth_data(
                        glob.glob('data/depth_files/' + mouse + '_' + tissue + '_' + 
                                  '*.dcs.region.mutpos.vcf_depth.txt')[0])
                    clone_freq = clone_count / median_depth
                    clone_type_freq = clone_type_count.div(median_depth)
    
                    temp_df = pd.DataFrame([mouse, group, tissue, clone_count,
                                            clone_type_count['AT_TA'][0],
                                            clone_type_count['AC_TG'][0],
                                            clone_type_count['AG_TC'][0],
                                            clone_type_count['GA_CT'][0],
                                            clone_type_count['GT_CA'][0],
                                            clone_type_count['GC_CG'][0],
                                            percent_clone, clone_freq,
                                            clone_type_freq['AT_TA'][0],
                                            clone_type_freq['AC_TG'][0],
                                            clone_type_freq['AG_TC'][0],
                                            clone_type_freq['GA_CT'][0],
                                            clone_type_freq['GT_CA'][0],
                                            clone_type_freq['GC_CG'][0]],
                                           index=["Mouse_ID", "Cohort", "Tissue",
                                                  "Clone_Count", "A>T/T>A_Count",
                                                  "A>C/T>G_Count", "A>G/T>C_Count",
                                                  "G>A/C>T_Count", "G>T/C>A_Count",
                                                  "G>C/C>G_Count", "Percent_Clone",
                                                  "Clone_Freq", "A>T/T>A_Freq",
                                                  "A>C/T>G_Freq", "A>G/T>C_Freq",
                                                  "G>A/C>T_Freq", "G>T/C>A_Freq",
                                                  "G>C/C>G_Freq"]).T
    
                    final_clone_data = pd.concat([final_clone_data, temp_df], 
                                                 axis=0)
                except:
                    pass
                
    final_clone_data.reset_index(drop=True, inplace=True)

    return final_clone_data


if __name__ == "__main__":

    # Import data
    if not os.path.isfile("data/imported_data/summary_data_wide.csv"):
        if not os.path.isdir("data/imported_data/"):
            os.mkdir("data/imported_data/")

        summary_data = summary_import('data/Mouse_aging_mtDNA_summary.csv')
        summary_data.to_csv("data/imported_data/summary_data_wide.csv")
    else:
        summary_data = pd.read_csv("data/imported_data/summary_data_wide.csv")

    if not os.path.isfile('data/imported_data/summary_data_tidy.csv'):
        summary_data_long = melt_summary(summary_data)
        summary_data_long.to_csv("data/imported_data/summary_data_tidy.csv")

    else:
        summary_data_long = pd.read_csv("data/imported_data/summary_data_tidy.csv")

    if not os.path.isfile("data/imported_data/mut_file_data.csv"):
        mut_data = mut_file_import()
        mut_data.to_csv("data/imported_data/mut_file_data.csv")
    else:
        mut_data = pd.read_csv("data/imported_data/mut_file_data.csv",
                               index_col=[0, 1])

    if not os.path.isfile("data/imported_data/summary_clone_data.csv"):
        final_clone_data = calc_clone_numbers(mut_data)
        final_clone_data.to_csv("data/imported_data/summary_clone_data.csv")
    else:
        final_clone_data = pd.read_csv("data/imported_data/summary_clone_data.csv")
