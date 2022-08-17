# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats

#def setup_figure():
def import_data(path1, path2):
    data1 = pd.read_csv(path1)
    data2 = pd.read_csv(path2)
    
    data1 = data1[["Age", "Treatment", "Tissue", "Total_SNV_Freq"]]
    data2 = data2[["Age", "Treatment", "Tissue", "Total_SNV_Freq"]]
    
    final_df = pd.concat([data1, data2])
    final_df.sort_values("Age", inplace=True)
    
    return final_df


def setup_figure():
    fig, ax = plt.subplots(nrows=3, dpi=600, figsize=(4,10), gridspec_kw={'hspace':0.5})

    return fig, ax    
    
  
def plot_subfig(data, tissue_list, title, ax):
    sub_data = data.query("Tissue in @tissue_list & Treatment=='NT'")
    
    color_vector = ["black"] * len(sub_data.query("Age==0.66")) + \
                   ["purple"] * len(sub_data.query("Age==4.5")) + \
                   ["black"] * len(sub_data.query("Age==10")) + \
                   ["purple"] * len(sub_data.query("Age==26"))
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(sub_data['Age'], 
                                                                   sub_data['Total_SNV_Freq'])
    
    plot = sns.regplot(x='Age', y='Total_SNV_Freq', data=sub_data, ax=ax, 
                       scatter_kws={"color":color_vector, "alpha":0.5}, 
                       line_kws={"color":"black", \
                                 'label':"$R^2={0:.2f}$\nP={1:.2e}".format(r_value**2, p_value)})
     
    plot.legend(ncol=1, frameon=False)    
    plot.set_title(title)
    
    plot.set_xlabel("Age (Months)")
    fig.canvas.draw()
    plt.setp(plot.get_yaxis().get_offset_text(), visible=False)
    order_of_mag = plot.get_yaxis().get_offset_text().get_text()[-2:]
    string = "Mutation Frequency ($\mathregular{10^{" + str(order_of_mag) + "}}$)"
    plot.set_ylabel(string)
    
    return plot

if __name__ == "__main__":
    data = import_data('data/misc_items/Arbeithuber_data_summary.csv', 'data/Mouse_aging_mtDNA_summary.csv')
    
    fig, ax = setup_figure()
    
    for index, tissues in enumerate([['M'], ['Br', 'C', 'Hi'], ['L', 'Li']]):
        title = ["Skeletal Muscle", "Brain", "Liver"][index]
        plot_subfig(data, tissues, title, ax[index])

    if not os.path.isdir("figures/"):
        os.mkdir("figures")

    fig.savefig('figures/Supp_Figure_5.png', dpi=600, facecolor='white', 
                bbox_inches='tight')
    fig.savefig('figures/Supp_Figure_5.pdf', dpi=600, facecolor='white', 
                bbox_inches='tight')