#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 15:58:36 2022

@author: scottrk
"""
print("Compiling and parsing raw data to tables for plotting. Compiled data are in 'data/imported data/'")
exec(open("compile_data.py").read())

print("Generating statistics and p-values. Per-figure statistics can be found in 'data/stats/'")
exec(open("compute_stats.py").read())

print("Generating Figure 1.")
exec(open("Figure_1.py").read())

print("Generating Figure 2.")
exec(open("Figure_2.py").read())

print("Generating Figure 3.")
exec(open("Figure_3.py").read())

print("Generating Figure 4")
exec(open("Figure_4.py").read())

print("Generating Figure 5")
exec(open("Figure_5.py").read())

print("Generating Supplemental Figure 1")
exec(open("Supp_Figure_1.py").read())

print("Generating Supplemental Figure 2")
exec(open("Supp_Figure_2.py").read())

print("Generating Supplemental Figure 3")
exec(open("Supp_Figure_3.py").read())

print("Generating Supplemental Figure 4")
exec(open("Supp_Figure_4.py").read())

print("Generating Supplemental Figure 8")
exec(open("Supp_Figure_7.py").read())

print("Generating Supplemental Figure 9")
exec(open("Supp_Figure_8.py").read())

print("Generating Supplemental Figure 10")
exec(open("Supp_Figure_9.py").read())

print("Generating Supplemental Figure 11")
exec(open("Supp_Figure_10.py").read())
