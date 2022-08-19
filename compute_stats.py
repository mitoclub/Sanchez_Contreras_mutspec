#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 09:36:16 2022

@author: scottrk
"""
import os
from collections import defaultdict
import scipy.stats as stats
import pandas as pd
from compile_data import summary_import, melt_summary, mut_file_import, \
                         calc_clone_numbers
from GlobalVars_ import tissue_type, mut_type, mut_type_conv, tissue_type_abbrev, \
                        R_lib_path
from HelperFuncs_ import statsmodel_summary_to_df
import numpy as np
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import rpy2.robjects as robjects


def Fig1A_D_stats(data, mut_class, output=None):
    
    p_value_dict = defaultdict(lambda:[])
    
    for tissue in tissue_type_abbrev[: -1]:
        young_data = data.query("Age=='Young' & Tissue==@tissue & Treatment=='NT'")
        old_data = data.query("Age=='Old' & Tissue==@tissue & Treatment=='NT'")
        
        F_stat, p_value = stats.ttest_ind(young_data[mut_class], old_data[mut_class], 
                                          equal_var=False)
        p_value_dict[tissue].append(p_value)
    
    if output != None:
        
        if os.path.isdir("data/stats/") == False:
            os.mkdir("data/stats/")
        
        pvalues = pd.DataFrame(p_value_dict).to_csv(output)
        
    return p_value_dict


def anova_df(df, cohort, mut_class):

    data = []

    for tissue in tissue_type_abbrev[: -1]:
        sub_data = df.query("Age==@cohort & Treatment=='NT' & Tissue==@tissue & Class==@mut_class")
        data.append(sub_data['Frequency'])

    f_stat, p_val = stats.f_oneway(data[0], data[1], data[2], data[3], data[4], 
                                   data[5], data[6], data[7])
    
    tukey_data = df.query("Age==@cohort & Treatment=='NT' & Class==@mut_class")
    tukey = pairwise_tukeyhsd(endog=tukey_data['Frequency'], 
                              groups=tukey_data['Tissue'])

    return p_val, tukey


def heatmap_stats(data, cohort, mut_class, output=None):

    tukey_p_val = defaultdict(lambda:[])

    p_val, anova = anova_df(data, cohort, mut_class)
    statsmodel_df = statsmodel_summary_to_df(anova)
    
    for tissue1 in tissue_type_abbrev[:-1]:
        
        for tissue2 in tissue_type_abbrev[:-1]:
            
            df = statsmodel_df.query("group1==@tissue1 | group2==@tissue1")
            x = df.query("group1==@tissue2 | group2==@tissue2")['p-adj']
            
            if tissue1 != tissue2:
                tukey_p_val[tissue1].append(float(x))
            else:
                tukey_p_val[tissue1].append(np.nan)
    
    pvals_df = pd.DataFrame(tukey_p_val)
    pvals_df.index = tissue_type[:-1]
    pvals_df = pvals_df.T
    
    if output != None:
        if os.path.isdir("data/stats/") == False:
            os.mkdir("data/stats/")
        pvals_df.to_csv(output)
        
    return pvals_df
  
    
def Fig2C_stats(lib_loc):
    
    r_source = robjects.r['source']
    
    if os.path.isdir('data/stats/') ==  False:
        os.mkdir('data/stats/')
    #output = open("data/stats/Figure_2_ratio_statistics.csv", 'w')
        
    r_source('fold_changes.R')
    
    #output.close()


def Fig3A_B_stats(data, column_name, output):
    
    clone_percent_pvals = {tissue: 0 for tissue in tissue_type[: -1]}

    for tissue in tissue_type_abbrev[: -1]:
        
        old = stats.ttest_ind(data.query("Tissue==@tissue & Cohort=='Old'")[column_name], 
                              data.query("Tissue==@tissue & Cohort=='Young'")[column_name], 
                              equal_var=False) 
        clone_percent_pvals[tissue] = old.pvalue
        pvalues = pd.DataFrame(clone_percent_pvals, columns=tissue_type[:-1], 
                               index=["Clone Percent P-values"])
        
        if os.path.isdir("data/stats/") == False:
            os.mkdir("data/stats/")
        pvalues.to_csv(output)
    
    return pvalues


def Fig5_stats(lib_loc):

    r_source = robjects.r('source')
    
    if os.path.isdir('data/stats/') ==  False:
        os.mkdir('data/stats/')
    #output = open("data/stats/Figure_5_Dunnett_statistics.csv", 'w')
        
    r_source('Dunnett_test.R')
    
    #output.close()


def Supp_Fig_7_stats(clone_freq_df, output):
    
    pvalue_list = []

    for tissue in ['K', 'L', 'Hi', 'C', 'M', 'B']:
        test = stats.ttest_ind(clone_freq_df.query("Tissue==@tissue & Cohort=='Old'")['Clone_Freq'], 
                               clone_freq_df.query("Tissue==@tissue & Cohort=='perf'")['Clone_Freq'], 
                               equal_var=False )
        pvalue_list.append(test.pvalue)

    pvalues = pd.DataFrame(pvalue_list, index=['K', 'L', 'Hi', 'C', 'M','B'], columns=['P value'])
    
    if not os.path.isdir("data/stats/"):
        os.mkdir("data/stats/")
    
    pvalues.to_csv(output)
    
    return pvalues


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
        final_clone_data = pd.read_csv("data/imported_data/summary_clone_data.csv",
                                       index_col=0)
    
    
    # Figure 1 statistics
    Fig1A_D_stats(summary_data, "Total_SNV_Freq", "data/stats/Figure_1A_statistics.csv")
    Fig1A_D_stats(summary_data, "Total_InDel_Freq", "data/stats/Figure_1D_statistics.csv")
    
    heatmap_stats(summary_data_long, "Young", 'Total_SNV_Freq', "data/stats/Figure_1A_Young_heatmap_stats.csv")
    heatmap_stats(summary_data_long, "Old", 'Total_SNV_Freq', "data/stats/Figure_1A_Old_heatmap_stats.csv")
    heatmap_stats(summary_data_long, "Young", 'Total_InDel_Freq', "data/stats/Figure_1D_Young_heatmap_stats.csv")
    heatmap_stats(summary_data_long, "Old", 'Total_InDel_Freq', "data/stats/Figure_1D_Old_heatmap_stats.csv")
    
    # Figure 2A-B statistics
    for age in ['Young', 'Old']:

        for i, mut_class in enumerate(mut_type):

            heatmap_stats(summary_data_long, age, mut_class, "data/stats/Figure_2_" + 
                          age + "_" + mut_type_conv[mut_class] + "_heatmap_stats.csv")

    # Figure 2C statistics
    Fig2C_stats(R_lib_path) # change package location as necessary
    
    # Figure 3 statistics
    Fig3A_B_stats(final_clone_data, "Percent_Clone", "data/stats/Figure_3A_statistics.csv")        
    Fig3A_B_stats(final_clone_data, "Clone_Freq", "data/stats/Figure_3B_statistics.csv")
    
    #Figure 5 statistics
    Fig5_stats(R_lib_path)
    
    #Supplemental Figure 8 statistics
    Supp_Fig_7_stats(final_clone_data, "data/stats/Supplemental_Figure_7_statistics.csv")