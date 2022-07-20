#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 13:25:12 2022

@author: scottrk
"""
cohort = ['Young', 'Old', 'Elam', 'NMN', 'perf']

old_mouse_id = ('JW24-R_NT', 'JW26-NP_NT', 'OM7-MS_NT', 'OM8-MS_NT', 'JW24-LF_NT', 'JW24-NP_NT')
young_mouse_id = ('YM13-MS_NT', 'YM14-MS_NT', 'YM3-MS_NT', 'YM5-MS_NT', 'YM7-MS_NT')
old_ss31_id = ('JW21-R_SS31', 'JW21-NP_SS31', 'JW55-RB_SS31', 'JW22-LF_SS31', 'JW22-NP_SS31')
old_nmn_id = ('JW27-LB_NMN', 'JW29-R_NMN', 'JW30-LF_NMN')
perfused_mouse_id = ('OM3-MS-P_perf', 'OM4-MS-P_perf', 'OM6-MS-P_perf')

#Tissues
tissue_type = ('K', 'L', 'EC', 'R', 'Hi', 'C', 'M', 'He', 'B')
tissue_type_long = ('Kidney', 'Liver',  'RPE/Chor.', 'Retina', 'Hippo.', 'Cerebel.', 'Sk. Muscle', 'Heart')

#Other
color_cycle = ['#ff9200', '#0433ff', '#ff84ff', '#76d5ff', '#932191', '#935200', '#008e00', '#ff2600']
mut_type = ['C>T/G>A', 'A>G/T>C', 'C>A/G>T', 'C>G/G>C','A>T/T>A', 'A>C/T>G']
mut_type_pretty = ['G→A/C→T', 'A→G/T→C', 'G→T/C→A', 'G→C/C→G', 'A→T/T→A', 'A→C/T→G']
mut_type_conv = {'C>T/G>A':'CT_GA', 'A>G/T>C':'AG_TC', 'C>A/G>T':'CA_GT', 'C>G/G>C':'CG_GC','A>T/T>A':'AT_TA','A>C/T>G':'AC_TG'}

R_lib_path =  '/home/scottrk/R/x86_64-pc-linux-gnu-library/4.1' # change package location as necessary

