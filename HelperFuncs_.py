#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 14:56:06 2022

@author: scottrk
"""
import numpy as np
import pandas as pd

def p_val_convert(df):
    
    df = np.where((df >= 0) & (df < 0.0001) ,4, df)
    df = np.where((df <= 1) & (df > 0.05), 0, df)
    df = np.where((0.0001 < df) & (df < 0.001), 3, df)
    df = np.where((0.001 < df) & (df < 0.01), 2, df)
    df = np.where((0.01 < df) & (df < 0.05), 1, df)
    
    return df


def statsmodel_summary_to_df(summary_table):
    
    summary_table_data = []
    
    for x in summary_table._results_table:
        summary_table_data.append(x.data)
        df = pd.DataFrame(summary_table_data[1:], columns=summary_table_data[0])
    
    return df