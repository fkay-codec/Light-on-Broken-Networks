"""
Here we are going to assess the normality of our edges across subjects.
"""

import os
import pandas as pd
import numpy as np
import scipy.stats as stats


def column_order_ismatch(fnirs_dict, fmri_dict):
    """
    
    Check if all dataframes in dict1 and dict2 have the same column order.
    
    Otherwise the correlation will be meaningless.
    """
    reference_columns = None
    # check fnirs dict
    for key, df in fnirs_dict.items():
        cols = df.columns.tolist()
        rows = df.index.tolist()
        if cols != rows:
            print(f"Warning: fnirs df {key} has non-symmetric labels (cols!=rows).")
        if reference_columns is None:
            reference_columns = cols
        elif cols != reference_columns:
            print(f"Column order mismatch in fnirs for key: {key}")
            return False

    # check fmri dict without mutating original frames
    for key, df in fmri_dict.items():
        # normalize the fmri column labels if they contain a prefix like "ROI1_S1D1"
        norm_cols = [c.split("_", 1)[-1] for c in df.columns.tolist()]
        norm_rows = [r.split("_", 1)[-1] for r in df.index.tolist()]
        if norm_cols != norm_rows:
            print(f"Warning: fmri df {key} has non-symmetric labels after normalization.")
        if norm_cols != reference_columns:
            print(f"Column order mismatch in fmri for key: {key}")
            print(f"Expected: {reference_columns}")
            print(f"Found:    {norm_cols}")
            return False

    return True

# Load correlation matrices for all subjects || Standard correlation matrices for testing
fnirs_folder = r"D:\Foivos\fNIRS_fMRI_Study\fNIRS_data\1_CorrelationMatrices\2_PartialCor"
fmri_folder = r"D:\Foivos\fNIRS_fMRI_Study\fMRI_data\3_CorrelationMatrices\2_PartialCor"

# Ask the user for input to compute for HbO or HbR
choice = input("Compute it for HbO or HbR? \nEnter 1 for HbO, 2 for HbR: ")

if choice == "1":
    modality = "HbO"
    print("You selected HbO.")
elif choice == "2":
    modality = "HbR"
    print("You selected HbR.")
else:
    print("Invalid input. Please enter 1 for HbO or 2 for HbR.")
    quit()

if modality == "HbO":
    excel_sheet_name = "HbO"
if modality == "HbR":
    excel_sheet_name = "HbR"

fnirs_all_data = {}
# load every subject's correlation matrix and place it in a dictionary
for file in os.listdir(fnirs_folder):
    if file.endswith("PartialCorM.xlsx"):
        subject_id = file.split("_")[0]
        fnirs_file = os.path.join(fnirs_folder, file)
        # load excel file
        fnirs_df = pd.read_excel(fnirs_file, index_col=0, sheet_name=excel_sheet_name)
        # clean column names
        fnirs_df.columns = [col.replace("_", "") for col in fnirs_df.columns]
        fnirs_df.index = [idx.replace("_", "") for idx in fnirs_df.index]
        fnirs_all_data[subject_id] = fnirs_df

fmri_all_data = {}
for file in os.listdir(fmri_folder):
    if file.endswith("PartialCorM.xlsx"):
        subject_id = file.split("_")[0]
        fmri_file = os.path.join(fmri_folder, file)
        # load excel file
        fmri_df = pd.read_excel(fmri_file, index_col=0)
        # no need to clean column names, now all column  names are compatible with each other
        fmri_all_data[subject_id] = fmri_df



if not column_order_ismatch(fnirs_all_data, fmri_all_data):
    print("Catastrophic error... Exiting...")
    quit()

##* I am not changing the below script, it was done for another purpose but i need the skeleton below, so it haas its use still...

# Take the fmri subject keys and use the first one to get the column names in order to create the dataframes with the results
fnirs_subjects = sorted(fnirs_all_data.keys())
fmri_subjects = sorted(fmri_all_data.keys())
fmri_single_subj = fmri_all_data[fmri_subjects[0]]
reference_col = fmri_single_subj.columns
reference_col = [c.replace("_", "/") for c in reference_col]

# Sanity check
if not column_order_ismatch(fnirs_all_data, fmri_all_data):
    print("Catastrophic error... Exiting...")
    quit()


#* not really needed but it is easier to leave it as is to perform the other things i want... normality check...
# We are only interested in the upper triangle of the resulting matrix, because the matrices are symmetric and the diagonal is always 1
# thus we will ignore the lower triangle and the diagonal
# in the results_df place from the fmri_df the columns
results_df = pd.DataFrame(index=reference_col, columns=reference_col)

######### Here we change stuff #########


txt_results = []
txt_results.append("="*100 + "\n")
txt_results.append("Assessing normality of edges across subjects using Shapiro-Wilk test\n")
txt_results.append("="*100 + "\n")
# create a variable to store how many times the null was accepted/rejected
fmri_shapirp_h0 = 0
fmri_shapirp_h1 = 0
fmri_tests_done = 0

fnirs_shapirp_h0 = 0
fnirs_shapirp_h1 = 0
fnirs_tests_done = 0



for i in range(len(results_df.index)): # rows
    for j in range(i+1, len(results_df.columns)): # columns
        fmri_edge = []
        fnirs_edge = []
        for fnirs_subj in fnirs_subjects:
            fnirs_edge.append(fnirs_all_data[fnirs_subj].iloc[i, j])
        for fmri_subj in fmri_subjects:
            fmri_edge.append(fmri_all_data[fmri_subj].iloc[i , j])

        r_clip = 0.999999  # to avoid inf values in Fisher transform
        fnirs_edge = np.clip(fnirs_edge, -r_clip, r_clip)
        fmri_edge = np.clip(fmri_edge, -r_clip, r_clip)
        # Fisher r-to-z transform
        fnirs_edge_z = np.arctanh(fnirs_edge)
        fmri_edge_z = np.arctanh(fmri_edge)
        # perform the shapiro-wilk test for normality
        sh_stat_fmri, p_val_fmri = stats.shapiro(fmri_edge_z)
        sh_stat_fnirs, p_val_fnirs = stats.shapiro(fnirs_edge_z)
        if p_val_fmri >= 0.05:
            fmri_shapirp_h0 += 1
        else:
            fmri_shapirp_h1 += 1

        if p_val_fnirs >= 0.05:
            fnirs_shapirp_h0 += 1
        else:
            fnirs_shapirp_h1 += 1
        fmri_tests_done += 1
        fnirs_tests_done += 1


txt_results.append(f"fMRI Shapiro-Wilk test results: \nH0 accepted (normality) {fmri_shapirp_h0} times\nH1 accepted (non-normality) {fmri_shapirp_h1} times\nTotal tests done: {fmri_tests_done}\n")
total_normal_fmri = (fmri_shapirp_h0/fmri_tests_done)*100
total_normal_fnirs = (fnirs_shapirp_h0/fnirs_tests_done)*100
txt_results.append(f"Total edges with normal distribution in fMRI: {total_normal_fmri:.2f}%\n")


txt_results.append("-"*50 + "\n")
txt_results.append(f"fNIRS ({modality}) Shapiro-Wilk test results: \nH0 accepted (normality) {fnirs_shapirp_h0} times\nH1 accepted (non-normality) {fnirs_shapirp_h1} times\nTotal tests done: {fnirs_tests_done}\n")
txt_results.append(f"Total edges with normal distribution in fNIRS ({modality}): {total_normal_fnirs:.2f}%\n")
txt_results.append("-"*50 + "\n\n")
txt_results.append("End of the results\n")
# save the results to a text file


# save the results
output_folder = r"D:\Foivos\fNIRS_fMRI_Study\1_CrossModality_EdgeComparison"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
output_file = os.path.join(output_folder, f"Edgewise_NormalityAssessment_ShapiroWilk_ParCor_{modality}.txt")
with open(output_file, "w") as f:
    f.writelines(txt_results)
quit()





######################################### End of Main Script #########################################


import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
import seaborn as sns
def plot_heatmap(df, title=str, save_path=None, community_labels=None):
    """
    Plot a heatmap for a 2D DataFrame (typically correlation matrices).
    
    Parameters:
    -----------
    df : pandas.DataFrame
        2D DataFrame to plot as heatmap
    title : str
        Title for the plot
    """
    ############# Plotting the p-values matrix ##################
    # Create a binary mask: 1 for p < 0.05 (significant), 0 for p >= 0.05 (non-significant)
    # True = 1, False = 0 
    binary_df = (df < 0.05).astype(int)

    # create mask for lower triangle
    mask = np.tril(np.ones_like(df, dtype=bool))
    # Define colors for the binary values
    colors = {0: "maroon", 1: "lightsteelblue"}  # 0 = non-significant, 1 = significant


    plt.figure(figsize=(10, 12))
    # Plot the heatmap
    sns.heatmap(binary_df, annot=False, square=True, linewidths=0.02, linecolor='lightgray',
                mask=mask, cbar=False, cmap=ListedColormap(colors.values()))
    

    # Add title and labels
    plt.title(title, fontsize=14)
    plt.xlabel("Channels", fontsize=10)
    plt.ylabel("Channels", fontsize=10)

    # Set all ticks to show every channel, positioned at the center of squares
    plt.xticks([i + 0.5 for i in range(len(df.columns))], df.columns, fontsize=5, rotation=90)
    plt.yticks([i + 0.5 for i in range(len(df.index))], df.index, fontsize=5)


    legend_labels = {0: "Similar (p ≥ 0.05)", 1: "Significantly Different (p < 0.05)"}
    # Add a custom legend
    legend_handles = [plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=color, markersize=10, label=legend_labels[label])
                      for label, color in colors.items()]

    # Add a custom legend
    legend_handles = [
        Line2D([0], [0], marker='s', color='w', markerfacecolor=color, markersize=10, label=legend_labels[label])
        for label, color in colors.items()
    ]

    # Add custom text to the legend
    legend_handles.append(Line2D([0], [0], color='w', label="FDR-corrected p-values\nT-tests on z-transformed Pearson r edges"))

    plt.legend(handles=legend_handles, bbox_to_anchor=(1.47, 0.5), loc='center right')
    
    

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")
     
    plt.show()




# Now we want to see which edges are significant/non significant, to do so we are building a custom plotting function
# Null p >= 0.05: mean connectivity strength for edge (i,j) is equal across groups == Success
# Alt p < 0.05: means differ. == Failure

# we are going to plot the p-value matrix only
plot_heatmap(p_values_fdr_df, title=f"fMRI vs fNIRS ({modality}) Edgewise T-Test P-Value Heatmap",
                save_path=os.path.join(output_folder, f"Edgewise_CrossModality_2samplesTtest_FDR_pvalues_StandCor_{modality}.png"))