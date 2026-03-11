"""
This script performs an edgewise t-test to compare connectivity strengths between two independent groups: 
fNIRS (HbO or HbR) and fMRI. 
It computes t-values and p-values for each edge (ROI pair) and applies FDR correction to control for multiple comparisons. 
The results include t-value matrices, FDR-corrected p-value matrices, and visualizations as heatmaps

Experimental Setup:
-	31 fNIRS subjects vs 31 fMRI subjects
-	Two independent groups (not the same people, not paired).
-	82 ROIs/channels per modality (assumed correspondence).
-	Individual-level: each subject -- one ROI×ROI correlation matrix. #! Fisher z-transform the correlation coefficients before any group-level analysis.
-	Group-level: average across subjects → one mean matrix per modality.

In this script we are going to perform an edgewise t-test between the two groups
    - What is edgewise comparison?
        - For each edge (i,j), we form two vectors:
            fMRI_edge = [r_subj1, r_subj2, ..., r_subjN]
            fNIRS_edge = [r_subj1, r_subj2, ..., r_subjN]
        - The t-test between fMRI_edge and fNIRS_edge quantifies inter-modality difference
        #! it is not possible to perform and edgewise correlation because the two groups are independent, and this would require paired data
    - We repeat this for all edges (i,j) to get a matrix of t-values and p-values.
    - We can then apply multiple comparison correction (e.g., FDR) to the p-values.

Interpretation:
    Null: mean connectivity strength for edge (i,j) is equal across groups.
    Alt: means differ.


INPUT:
    - fNIRS correlation matrices for all subjects || Standard correlation matrices
    - fMRI correlation matrices for all subjects || Standard correlation matrices
OUTPUT:
    - Matrix of t-values
    - Matrix of p-values
    - FDR-corrected p-values
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


# We are only interested in the upper triangle of the resulting matrix, because the matrices are symmetric and the diagonal is always 1
# thus we will ignore the lower triangle and the diagonal

# in the results_df place from the fmri_df the columns
results_df = pd.DataFrame(index=reference_col, columns=reference_col)
# There is no need to populate the diagonal with 1s... Instead with we populate by NaN
np.fill_diagonal(results_df.values, np.nan)
# same for the p-values
p_results_df = pd.DataFrame(index=reference_col, columns=reference_col)
np.fill_diagonal(p_results_df.values, np.nan)

for i in range(len(results_df.index)): # rows
    for j in range(i+1, len(results_df.columns)): # columns
        # get the labels for the current edge
        # row_label = results_df.index[i]
        # col_label = results_df.columns[j]
        # print(f"Computing correlation for edge ({row_label}, {col_label}) [{i},{j}]")
        # quit()
        fmri_edge = []
        fnirs_edge = []
        for fnirs_subj in fnirs_subjects:
            fnirs_edge.append(fnirs_all_data[fnirs_subj].iloc[i, j])
            # print(f"appended row label {fnirs_all_data[fnirs_subj].index[i]} col label {fnirs_all_data[fnirs_subj].columns[j]} value {fnirs_all_data[fnirs_subj].iloc[i, j]} from subject {fnirs_subj}")
        # quit()
        for fmri_subj in fmri_subjects:
            fmri_edge.append(fmri_all_data[fmri_subj].iloc[i , j])
            # print(f"appended row label {fmri_all_data[fmri_subj].index[i]} col label {fmri_all_data[fmri_subj].columns[j]} value {fmri_all_data[fmri_subj].iloc[i, j]} from subject {fmri_subj}")
        r_clip = 0.999999  # to avoid inf values in Fisher transform
        fnirs_edge = np.clip(fnirs_edge, -r_clip, r_clip)
        fmri_edge = np.clip(fmri_edge, -r_clip, r_clip)
        # Fisher r-to-z transform
        fnirs_edge_z = np.arctanh(fnirs_edge)
        fmri_edge_z = np.arctanh(fmri_edge)
        # perform the two samples indepdendent t-test
        t_stat, p_val = stats.ttest_ind(fmri_edge_z, fnirs_edge_z, equal_var=False)  # Welch's t-test
        print(f"Edge ({i},{j}) t-stat: {t_stat}, p-value: {p_val}")
        quit()
        results_df.iloc[i, j] = t_stat
        results_df.iloc[j, i] = np.nan  # We dont care about the lower triangle, doesnt add anything new
        p_results_df.iloc[i, j] = p_val
        p_results_df.iloc[j, i] = np.nan

# FDR correction across all edges
upper_indices = np.triu_indices(len(results_df.columns), k=1)  # get upper triangle indices, excluding diagonal
p_values_upper = p_results_df.values[upper_indices]

p_values_upper = np.array(p_values_upper, dtype=float)
# apply FDR correction
p_values_fdr = stats.false_discovery_control(p_values_upper, method='bh')

# Fill it back to the matrix
p_values_fdr_df = p_results_df.copy()
p_values_fdr_df.values[upper_indices] = p_values_fdr

# Set lower triangle AND diagonal to NaN
p_values_fdr_df.values[np.tril_indices(len(p_values_fdr_df), k=0)] = np.nan

# save the results
output_folder = r"D:\Foivos\fNIRS_fMRI_Study\1_CrossModality_EdgeComparison"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
with pd.ExcelWriter(os.path.join(output_folder, f"Edgewise_CrossModality_2samplesTtest_ParCor_{modality}.xlsx")) as writer:
    results_df.to_excel(writer, sheet_name='Edgewise_tvalues')
    p_values_fdr_df.to_excel(writer, sheet_name='FDR_pvalues')





######################################### End of Main Script #########################################

quit()
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
                save_path=os.path.join(output_folder, f"Edgewise_CrossModality_2samplesTtest_FDR_pvalues_ParCor_{modality}.png"))