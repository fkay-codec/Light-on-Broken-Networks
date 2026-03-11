import os
import pandas as pd
import matplotlib.pyplot as plt
import random
import numpy as np
import seaborn as sns
from bct import density_dir



"""
lets create a histogram/bar plot for data to see the distribution of the matrix values

Option A:
X-Axis: edge value bins
Y-Axis: counts
==> one for fMRI one for fNIRS HbO one for fNIRS HbR in the same plot 


"""

def make_histogram(fmri_excel_path, fnirs_excel_path, out_folder_figures):
        
    fmri_df = pd.read_excel(fmri_excel_path, index_col=0)

    # make the diagonal and lower triangle nan
    # Mask lower triangle and diagonal
    fmri_df.values[np.tril_indices_from(fmri_df)] = np.nan
    edge_values_fmri = fmri_df.values.flatten()
    edge_values_fmri = edge_values_fmri[~np.isnan(edge_values_fmri)]

    # do the same for fNIRS HbO just for now
    fnirs_df = pd.read_excel(fnirs_excel_path, index_col=0)
    fnirs_df.values[np.tril_indices_from(fnirs_df)] = np.nan
    edge_values_fnirs = fnirs_df.values.flatten()
    edge_values_fnirs = edge_values_fnirs[~np.isnan(edge_values_fnirs)]

    fnirs_hbr_df = pd.read_excel(fnirs_excel_path, index_col=0, sheet_name='HbR_Average')
    fnirs_hbr_df.values[np.tril_indices_from(fnirs_hbr_df)] = np.nan
    edge_values_fnirs_hbr = fnirs_hbr_df.values.flatten()
    edge_values_fnirs_hbr = edge_values_fnirs_hbr[~np.isnan(edge_values_fnirs_hbr)]


    plt.figure(figsize=(10, 6))
    # lets do the same plot with sns
    sns.histplot(data=edge_values_fmri, bins=250, color='tab:green', alpha=0.6, label='fMRI', kde=True)#, element='step' )
    sns.histplot(data=edge_values_fnirs, bins=250, color='tab:red', alpha=0.5, label='fNIRS HbO', kde=True)#, element='step')
    sns.histplot(data=edge_values_fnirs_hbr, bins=250, color='tab:blue', alpha=0.4, label='fNIRS HbR', kde=True)#, element='step')
    
    if 'Partial' in os.path.basename(fmri_excel_path):
        type_of_plot = 'Partial'
    elif 'Standard' in os.path.basename(fmri_excel_path):
        type_of_plot = 'Bivariate'
    title = f"Edge Weight Histogram\n ({type_of_plot} Correlations)"
    
    plt.title(title, fontsize=15, fontweight='bold')
    plt.xlabel('Edge Value (Fisher z)')
    plt.ylabel('Count')
    plt.legend()

    # plt.tight_layout()

    # save figure
    figure_name = f"Histogram_Average_{type_of_plot.replace(' ', '_')}.png"
    plt.savefig(os.path.join(out_folder_figures, figure_name), dpi=300)
    # plt.show()



out_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\2_CrossModalityGroupAnalysis"
out_folder_figures = os.path.join(out_folder, 'GroupMatrices_Figures')
if not os.path.exists(out_folder_figures):
    os.makedirs(out_folder_figures)

fmri_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\4_AverageCorM"
fmri_excel_paths = [f for f in os.listdir(fmri_folder) if f.endswith('.xlsx')]

fnirs_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\2_AverageCorM"
fnirs_excel_paths = [f for f in os.listdir(fnirs_folder) if f.endswith('.xlsx')]

for fmri_file in fmri_excel_paths:
    for fnirs_file in fnirs_excel_paths:
        if fmri_file == fnirs_file:
            print(f"Matched files: {fmri_file} and {fnirs_file}")
            # visualize these two files
            make_histogram(fmri_excel_path=os.path.join(fmri_folder, fmri_file), 
                           fnirs_excel_path=os.path.join(fnirs_folder, fnirs_file), 
                           out_folder_figures=out_folder_figures)
            
