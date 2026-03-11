"""
In this script for each file found in the input folder, the average correlation matrix will be computed
"""


import os
import pandas as pd
import numpy as np
def fisher_z_transform(corr_matrix):
    """Apply Fisher z-transformation to a correlation matrix"""
    clip = 0.9999999999  # To avoid infinities
    corr_matrix = corr_matrix.clip(-clip, clip)  # Avoid infinities
    z_matrix = np.arctanh(corr_matrix)
    return z_matrix

input_master_folder = r"D:\Foivos\fNIRS_fMRI_Study\fMRI_data\3_CorrelationMatrices"
output_master_folder = r"D:\Foivos\fNIRS_fMRI_Study\fMRI_data\4_AverageCorM"
if not os.path.exists(output_master_folder):
    os.makedirs(output_master_folder)

for subfolder in os.listdir(input_master_folder):
    channel_names = None
    matrices = []
    matrices_3d=None
    avg_matrix=None
    avg_df=None

    input_folder = os.path.join(input_master_folder, subfolder)
    for file in os.listdir(input_folder):
        if not file.endswith("CorM.xlsx"): # skip the non correlation matrix files
            continue
        # load the excel file to a dataframe
        input_file = os.path.join(input_folder, file)
        print(f"Processing file: {input_file}")
        data_df = pd.read_excel(input_file, index_col=0) 

        if channel_names is None:
            channel_names = list(data_df.columns)

        # Apply Fisher z-transformation
        data_df = fisher_z_transform(data_df)
        
        # append the matrix to the lsit
        matrices.append(data_df.values)
    # convert the list to a 3D array
    matrices_3d = np.array(matrices)
    # average across subjects the 3d array
    avg_matrix = np.mean(matrices_3d, axis=0)
    avg_df = pd.DataFrame(avg_matrix, columns=channel_names, index=channel_names)

    subf_name = subfolder.split("_")[-1]
    # Save the average correlation matrix to an excel
    out_avg_file = os.path.join(output_master_folder, f"Group_Average_{subf_name}_Ztransf.xlsx")
    with pd.ExcelWriter(out_avg_file) as writer:
        avg_df.to_excel(writer, sheet_name='AverageCorM')
    print(f"Average correlation matrix saved to: {out_avg_file}")

