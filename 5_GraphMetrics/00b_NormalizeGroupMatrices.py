"""
To compute the graph metrics we are going to normalize all the adjacency matrices to the same reference

"""

import os
import pandas as pd
import numpy as np

def normalize_fmri_df(fmri_path, global_max):
    df = pd.read_excel(fmri_path, index_col=0)
    matrix = df.to_numpy()
    # remove the diagonal and fill the global max there instead
    np.fill_diagonal(matrix, global_max)

    # normalize the matrix
    normalized_matrix = matrix / global_max
    
    # create a new dataframe
    normalized_df = pd.DataFrame(normalized_matrix, index=df.index, columns=df.columns)
    return normalized_df

def find_max_fmri_value(fmri_path):
    df = pd.read_excel(fmri_path, index_col=0)
    matrix = df.to_numpy()
    # remove the diagonal
    np.fill_diagonal(matrix, 0)

    # get the max absolute value
    max_value = np.max(np.abs(matrix))
    return max_value

def find_max_fnirs_value(fnirs_path):
    sheets = pd.ExcelFile(fnirs_path).sheet_names
    df_oxy = pd.read_excel(fnirs_path, index_col=0, sheet_name=sheets[0])
    matrix_oxy = df_oxy.to_numpy()
    # remove the diagonal
    np.fill_diagonal(matrix_oxy, 0)

    # get the max absolute value
    max_value_oxy = np.max(np.abs(matrix_oxy))

    df_deoxy = pd.read_excel(fnirs_path, index_col=0, sheet_name=sheets[1])
    matrix_deoxy = df_deoxy.to_numpy()
    # remove the diagonal
    np.fill_diagonal(matrix_deoxy, 0)
    max_value_deoxy = np.max(np.abs(matrix_deoxy))
    return max_value_oxy, max_value_deoxy

def normalize_fnirs_df(fnirs_path, global_max):
    sheets = pd.ExcelFile(fnirs_path).sheet_names
    df_oxy = pd.read_excel(fnirs_path, index_col=0, sheet_name=sheets[0])
    matrix_oxy = df_oxy.to_numpy()
    # remove the diagonal and fill the global max there instead
    np.fill_diagonal(matrix_oxy, global_max)

    # normalize the matrix
    normalized_matrix_oxy = matrix_oxy / global_max
    
    # create a new dataframe
    normalized_oxy_df = pd.DataFrame(normalized_matrix_oxy, index=df_oxy.index, columns=df_oxy.columns)

    df_deoxy = pd.read_excel(fnirs_path, index_col=0, sheet_name=sheets[1])
    matrix_deoxy = df_deoxy.to_numpy()
    # remove the diagonal and fill the global max there instead
    np.fill_diagonal(matrix_deoxy, global_max)

    # normalize the matrix
    normalized_matrix_deoxy = matrix_deoxy / global_max
    
    # create a new dataframe
    normalized_deoxy_df = pd.DataFrame(normalized_matrix_deoxy, index=df_deoxy.index, columns=df_deoxy.columns)

    return normalized_oxy_df, normalized_deoxy_df
# First we are going to do it on Group Average matrices.

fmri_folder_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\4_AverageCorM"
fnirs_folder_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\2_AverageCorM"

load_fmri_files=[f for f in os.listdir(fmri_folder_path) if f.endswith(".xlsx")]
load_fnirs_files=[f for f in os.listdir(fnirs_folder_path) if f.endswith(".xlsx")]


for fmri_file in load_fmri_files:
    fmri_matrix_path = os.path.join(fmri_folder_path, fmri_file)
    if "Partial" in fmri_file:
        fmri_max_abs_par = find_max_fmri_value(fmri_matrix_path)
    if "Standard" in fmri_file:
        fmri_max_abs_biv = find_max_fmri_value(fmri_matrix_path)
print(f"fMRI Bivariate Max Abs Value: {fmri_max_abs_biv}")
print(f"fMRI Partial Max Abs Value: {fmri_max_abs_par}")

for fnirs_file in load_fnirs_files:
    fnirs_matrix_path = os.path.join(fnirs_folder_path, fnirs_file)
    if "Partial" in fnirs_file:
        fnirs_max_abs_par_oxy, fnirs_max_abs_par_deoxy = find_max_fnirs_value(fnirs_matrix_path)
    if "Standard" in fnirs_file:
        fnirs_max_abs_biv_oxy, fnirs_max_abs_biv_deoxy = find_max_fnirs_value(fnirs_matrix_path)
print(f"fNIRS Bivariate Max Abs Value HbO: {fnirs_max_abs_biv_oxy}")
print(f"fNIRS Bivariate Max Abs Value HbR: {fnirs_max_abs_biv_deoxy}")
print(f"fNIRS Partial Max Abs Value HbO: {fnirs_max_abs_par_oxy}")
print(f"fNIRS Partial Max Abs Value HbR: {fnirs_max_abs_par_deoxy}")

# find the global maximum for partial and bivariate
global_max_bivariate = max(fnirs_max_abs_biv_oxy, fnirs_max_abs_biv_deoxy, fmri_max_abs_biv)
global_max_partial = max(fnirs_max_abs_par_oxy, fnirs_max_abs_par_deoxy, fmri_max_abs_par)
print(f"Global Max Bivariate: {global_max_bivariate}")
print(f"Global Max Partial: {global_max_partial}")

out_for_fmri = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\5_AverCorM_NormForGraph"
if not os.path.exists(out_for_fmri):
    os.makedirs(out_for_fmri)


# now normalize all the matrices to their corresponding global max
for fmri_file in load_fmri_files:
    fmri_matrix_path = os.path.join(fmri_folder_path, fmri_file)
    normalized_df = pd.DataFrame()
    if "Partial" in fmri_file:
        normalized_df = normalize_fmri_df(fmri_matrix_path, global_max_partial)
    if "Standard" in fmri_file:
        normalized_df = normalize_fmri_df(fmri_matrix_path, global_max_bivariate)
    out_path = os.path.join(out_for_fmri, fmri_file.split('.')[0] + "_Normalized.xlsx")
    with pd.ExcelWriter(out_path) as writer:
        normalized_df.to_excel(writer, index=True)
# Save the normalized fMRI Dataframes



# the same for fNIRS
out_for_fnirs = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\3_AverCorM_NormForGraph"
if not os.path.exists(out_for_fnirs):
    os.makedirs(out_for_fnirs)

for fnirs_file in load_fnirs_files:
    fnirs_matrix_path = os.path.join(fnirs_folder_path, fnirs_file)
    normalized_oxy_df = pd.DataFrame()
    normalized_deoxy_df = pd.DataFrame()
    if "Partial" in fnirs_file:
        normalized_oxy_df, normalized_deoxy_df = normalize_fnirs_df(fnirs_matrix_path, global_max_partial)
    if "Standard" in fnirs_file:
        normalized_oxy_df, normalized_deoxy_df = normalize_fnirs_df(fnirs_matrix_path, global_max_bivariate)
    out_path = os.path.join(out_for_fnirs, fnirs_file.split('.')[0] + "_Normalized.xlsx")
    with pd.ExcelWriter(out_path) as writer:
        normalized_oxy_df.to_excel(writer, index=True, sheet_name="HbO")
        normalized_deoxy_df.to_excel(writer, index=True, sheet_name="HbR")
print("Normalization Complete!")
