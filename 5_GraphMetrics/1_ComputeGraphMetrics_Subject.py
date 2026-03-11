
'''
Script: 1_ComputeGrMetSubj.py

Description:
    This script computes node-level (nodal) graph metrics for fNIRS and fMRI correlation matrices.
    Global network metrics are NOT computed in this script (computed in group-level script instead).
    Results are saved in Excel files for further analysis.

Inputs:
    - fNIRS Correlation Matrices:
        Folder: "fnirs_folder" 
        File Format: Excel files (.xlsx) with sheets for HbO and HbR correlation matrices.
    
    - fMRI Correlation Matrices:
        Folder: "fmri_folder" 
        Structure: Subfolders for "1_StandardCor" and "2_PartialCor" correlation matrices.
        File Format: Excel files (.xlsx) with single sheet per subject.
    
    NOTE: Input data are Fisher Z-transformed correlation matrices that were then normalized 
    to the range [0, 1] by dividing by the global maximum absolute value across ALL subjects 
    and modalities. This enables direct comparison of graph metrics between fMRI and fNIRS.

Outputs:
    - Processed Graph Metrics:
        Folder: "fnirs_output_folder" and "fmri_output_folder"
        Structure: Subfolders matching the input folder structure (1_StandardCor/, 2_PartialCor/)
        File Format: Excel files (.xlsx) with suffix "_GraphMetrics.xlsx"
        - fNIRS: One sheet per hemoglobin type (HbO_GraphMetrics, HbR_GraphMetrics)
        - fMRI: Single sheet named "GraphMetrics"

Metrics Computed (Node-Level Only):
    1. Nodal Strength Variants:
        - "Net Nodal Strength": Total weighted degree (sum of all connections)
        - "Positive Nodal Strength": Sum of positive connections only
        - "Negative Nodal Strength": Sum of negative connections only
        - "Absolute Nodal Strength": Sum of positive + absolute value of negative
        - "Normalized Nodal Strength": Rubinov & Sporns (2011) formula accounting for pos/neg balance
    
    2. Local Assortativity:
        - "Positive Assortativity": Local assortativity using positive weights
        - "Negative Assortativity": Local assortativity using negative weights
        - Computed using brainconn library (compatible with signed weights)
    
    3. Local Efficiency:
        - "Local Efficiency": Node-level efficiency (BCT algorithm)
        - Uses absolute values of weights (negatives treated as positive)
        - BCT requires positive weights in range [0, 1]
    
    4. Nodal Path Length:
        - "Nodal Path Length [1-|z|]": Average shortest path per node
        - Distance metric: 1 - |z| (treats positive and negative correlations equally)
        - Excludes infinite distances (disconnected nodes) from average
    
    5. Betweenness Centrality:
        - "Betweenness Centrality [1-|z|]": Node importance based on shortest paths
        - Uses same distance metric as path length: 1 - |z|

Libraries Used:
    - BCT (Brain Connectivity Toolbox): Path length, betweenness, local efficiency, nodal strength
    - brainconn: Local assortativity with signed weights

Additional Notes:
    - Script preserves signed weights (positive and negative) where applicable
    - Distance metric 1-|z| treats strong negative correlations as equally "close" to strong positives
    - Diagonal values are set to 0 before all computations
'''

import warnings


import networkx
import numpy as np
import pandas as pd
import bct.algorithms as bct
import os
import brainconn.core as bc

# Suppress RuntimeWarnings globally
warnings.filterwarnings("ignore", category=RuntimeWarning)


def nodal_path_length(D):
    nodal_path_length = np.zeros(D.shape[0])

    # get the average distance while ignoring inf values (disconnected nodes), and self
    for i in range(D.shape[0]):
        mask = (np.arange(D.shape[0]) != i) & (D[i, :] != np.inf)
        nodal_path_length[i] = np.mean(D[i, mask])
    return nodal_path_length

def nodal_strength(df):
    df=df.copy()
    array = df.to_numpy()
    np.fill_diagonal(array, 0)
    df = pd.DataFrame(array, index=df.index, columns=df.columns)
    node_str = df.sum(axis=1) # all columns for each row
    return node_str

def compute_graph_metrics(corr_df):
    """
    Compute node-level graph metrics from a normalized correlation matrix.
    
    Input: 
        corr_df: DataFrame containing Fisher Z-transformed correlation matrix that has been 
                normalized to range [0, 1] by dividing by global maximum across all subjects/modalities.
                Can be from individual subject (not group average).
    
    Output: 
        DataFrame with node-level metrics (one row per node/channel):
        - Net/Positive/Negative/Absolute/Normalized Nodal Strength
        - Positive/Negative Local Assortativity
        - Local Efficiency
        - Nodal Path Length using distance metric 1-|z|
        - Betweenness Centrality using distance metric 1-|z|
    """

    corr_array = corr_df.to_numpy().copy()
    # Remove the diagonal and fill with 0
    np.fill_diagonal(corr_array, 0)


    # Compute nodal strength variants
    # Data is Fisher Z-transformed and globally normalized, enabling cross-modal comparison
    # Global normalization preserves mean differences between groups (unlike within-matrix z-scoring)

    # net node str
    node_str = bct.strengths_und(corr_array)
    # compute the nodal str for pos and neg
    node_str_pos, node_str_neg, _, _ = bct.strengths_und_sign(corr_array)
    
    
    # Compute normalized nodal strength based on Rubinov and Sporns (2011)
    # Formula accounts for balance between positive and negative connections
    # Dividing by 81 normalizes by (N-1) where N=82 nodes
    norm_node_str = node_str_pos/81 - (node_str_neg/81)*(node_str_pos/(node_str_pos +node_str_neg))

    # Compute the absolute nodal strength
    node_str_abs = node_str_pos + np.abs(node_str_neg)

    print("             Computing assortativity...")
    # Compute local assortativity using brainconn library (requires diagonal = 0)
    # This function handles signed weights (positive and negative correlations)
    pos_assort, neg_assort = bc.local_assortativity_wu_sign(corr_array)


    # Compute shortest paths and betweenness centrality
    # Requires converting correlation matrix to distance matrix
    # Our matrices are in range [0, 1] after Fisher Z-transform and global normalization
    
    z = corr_df.to_numpy().copy()

    # Create distance matrix using: L = 1 - |z|
    # This treats strong negative correlations as equally "close" to strong positive correlations
    # (e.g., both r=0.9 and r=-0.9 would have distance 0.1)
    tmp = z.copy()
    tmp = np.abs(tmp)
    L1 = 1 - tmp


    print("             Creating distance matrix...")
    # Compute all-pairs shortest path distances using BCT
    D1, _ = bct.distance_wei(L1)

    print("             Computing path lengths...")
    # Compute average shortest path length per node (excluding self and infinite distances)
    nodal_path_length1 = nodal_path_length(D1)



    print("             Computing betweenness centrality...")
    # Compute betweenness centrality per node using distance metric 1-|z|
    # Measures importance of each node in shortest paths between other nodes
    nodal_betweeness1 = bct.betweenness_wei(L1)

    print("             Computing local efficiency...")
    # Local efficiency using BCT algorithm
    # BCT requires positive weights in range [0, 1], so we use absolute values
    # This treats negative correlations as positive with same magnitude
    tmp = z.copy()
    tmp = np.abs(tmp)
    local_eff = bct.efficiency_wei(tmp, local=True)


    print("            Preparing DataFrame...")
    # convert to pandas series for easier handling and saving
    # node str
    node_str = pd.Series(node_str, index=corr_df.index)
    node_str_pos = pd.Series(node_str_pos, index=corr_df.index)
    node_str_neg = pd.Series(node_str_neg, index=corr_df.index)
    node_str_abs = pd.Series(node_str_abs, index=corr_df.index)
    norm_node_str=pd.Series(norm_node_str, index=corr_df.index)
    # assortativity
    pos_assort = pd.Series(pos_assort, index=corr_df.index)
    neg_assort = pd.Series(neg_assort, index=corr_df.index)

    # path length
    nodal_path_length1 = pd.Series(nodal_path_length1, index=corr_df.index)

    # betweeness centrality
    nodal_betweeness1 = pd.Series(nodal_betweeness1, index=corr_df.index)
    
    # local efficiency
    local_eff = pd.Series(local_eff, index=corr_df.index)

    dataframe = pd.DataFrame({
        "Net Nodal Strength": node_str,
        "Positive Nodal Strength": node_str_pos,
        "Negative Nodal Strength": node_str_neg,
        "Absolute Nodal Strength": node_str_abs,
        "Normalized Nodal Strength": norm_node_str,
        "Positive Assortativity": pos_assort,
        "Negative Assortativity": neg_assort,
        "Local Efficiency": local_eff,
        "Nodal Path Length [1-|z|]": nodal_path_length1,
        "Betweenness Centrality [1-|z|]": nodal_betweeness1,
    }, index=corr_df.index)    
    # Write the DataFrame to the Excel file
    return dataframe


# Load the folders with the data
fnirs_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\1b_CorrelationMatrices_NormForGraph"
fnirs_output_folder = os.path.dirname(fnirs_folder) + r"\4_GraphMetricsPerSubject"
if not os.path.exists(fnirs_output_folder):
    os.makedirs(fnirs_output_folder)
fmri_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\3b_CorrelationMatrices_NormForGraph"
fmri_output_folder = os.path.dirname(fmri_folder) + r"\6_GraphMetricsPerSubject"
if not os.path.exists(fmri_output_folder):
    os.makedirs(fmri_output_folder)

print("Processing fNIRS data...")
# Process both StandardCor and PartialCor folders

for folder in os.listdir(fnirs_folder):
    if "1_Standard" not in folder and "2_Par" not in folder:
        continue

    print(f"Processing folder: {folder}")
    fnirs_subfolder = os.path.join(fnirs_folder, folder)
    # lets load the excel files for each heme group
    for file in os.listdir(fnirs_subfolder):
        if not file.endswith("xlsx"):
            continue
        print(f"    Processing file: {file}")
        fnirs_file = os.path.join(fnirs_subfolder, file)
        sheet_names = pd.ExcelFile(fnirs_file).sheet_names

        output_folder = os.path.join(fnirs_output_folder, folder)
        output_file = os.path.join(output_folder, f"{file.split('.')[0]}_GraphMetrics.xlsx")

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            print(sheet_names)
            for sheet in sheet_names:
                print(f"        Processing sheet: {sheet}")
                fnirs_df = pd.read_excel(fnirs_file, sheet_name=sheet, index_col=0)

                fnirs_array = fnirs_df.to_numpy().copy()
                print("         Computing graph metrics...")
                GraphMetrics_df = compute_graph_metrics(fnirs_df)

                # Save the graph metrics to the Excel file
                # Write the DataFrame to the Excel file
                GraphMetrics_df.to_excel(writer, sheet_name=f"{sheet}_GraphMetrics")
print("Finished processing fNIRS data.\n\n")


print("Now processing fMRI data...")
# the same for fMRI data
for folder in os.listdir(fmri_folder):
    if "1_Standard" not in folder and "2_Par" not in folder:
        continue
    print(f"Processing folder: {folder}")
    fmri_subfolder = os.path.join(fmri_folder, folder)
    # Load the Excel files for each group
    for file in os.listdir(fmri_subfolder):
        if not file.endswith("xlsx"):
            continue
        print(f"        Processing file: {file}")
        fmri_file = os.path.join(fmri_subfolder, file)

        output_folder = os.path.join(fmri_output_folder, folder)
        output_file = os.path.join(output_folder, f"{file.split('.')[0]}_GraphMetrics.xlsx")

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        fmri_df = pd.read_excel(fmri_file, index_col=0)
        # No need for Fisher Z-transform here - already done in normalization script (00a)
        
        fmri_array = fmri_df.to_numpy().copy()
        print("         Computing graph metrics...")
        GraphMetrics_df = compute_graph_metrics(fmri_df)

        # save the graph metric to the excel file
        
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            GraphMetrics_df.to_excel(writer, sheet_name="GraphMetrics")


print(f"Finished Processing")