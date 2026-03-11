"""adapted from partial correlation folder"""
import os
import pandas as pd
import utils as ut
import numpy as np
import scipy.stats as stats

def fisher_z_transform(corr_matrix):
    """Apply Fisher z-transformation to a correlation matrix"""
    clip = 0.9999999999  # To avoid infinities
    corr_matrix = corr_matrix.clip(-clip, clip)  # Avoid infinities
    z_matrix = np.arctanh(corr_matrix)
    return z_matrix



input_folder = r"D:\Foivos\fNIRS_fMRI_Study\fNIRS_data\1_CorrelationMatrices\1_StandardCor"
output_folder = r"D:\Foivos\fNIRS_fMRI_Study\fNIRS_data\2_AverageCorM"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Step 1: Load all matrices into 3D arrays
hbo_matrices = []
hbr_matrices = []
channel_names = None



for file in os.listdir(input_folder):
    if not file.endswith(".xlsx"): # skip the non excel files
        continue

    input_file = os.path.join(input_folder, file)
    print(f"Processing file: {input_file}")
    
    # Load the HbO sheet from the Excel file
    hbo_df = pd.read_excel(input_file, sheet_name='HbO', index_col=0)  # Use first column as index
    hbr_df = pd.read_excel(input_file, sheet_name='HbR', index_col=0)  # Use first column as index

    # Store channel names from the first file
    if channel_names is None:
        channel_names = list(hbo_df.columns)

    # z-transform the correlation matrices
    hbo_df = fisher_z_transform(hbo_df)
    hbr_df = fisher_z_transform(hbr_df)

    # Append matrices to the lists
    hbo_matrices.append(hbo_df.values)
    hbr_matrices.append(hbr_df.values)

# Convert lists to 3D arrays: (n_subjects, n_channels, n_channels)
hbo_3d = np.array(hbo_matrices)
hbr_3d = np.array(hbr_matrices)

# Vectorized average computation - much simpler than correlation!
hbo_avg_matrix = np.mean(hbo_3d, axis=0)  # Average across subjects (axis=0)
hbr_avg_matrix = np.mean(hbr_3d, axis=0)  # Average across subjects (axis=0)

print(f"Average computation completed!")
print(f"HbO average matrix shape: {hbo_avg_matrix.shape}")
print(f"HbR average matrix shape: {hbr_avg_matrix.shape}")

# Step 4b: Store average values in new matrix maintaining same structure
hbo_avg_df = pd.DataFrame(hbo_avg_matrix, columns=channel_names, index=channel_names)
hbr_avg_df = pd.DataFrame(hbr_avg_matrix, columns=channel_names, index=channel_names)



output_file_avg = os.path.join(output_folder, "Group_Average_StandardCor_Ztransf.xlsx")

with pd.ExcelWriter(output_file_avg) as writer:
    hbo_avg_df.to_excel(writer, sheet_name='HbO_Average')
    hbr_avg_df.to_excel(writer, sheet_name='HbR_Average')

