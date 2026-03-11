"""
To compute the graph metrics we are going to normalize all the adjacency matrices to the same reference


Across subjects what we will do is the following:

for Partial:
first fisher r to z all the subject matrices
then

search for max across fNIRS HbO
search for max across fNIRS HbR
search for max across fMRI
get the global max

then normalize all partial correlation matrices to that global max
"""

import os
import pandas as pd
import numpy as np
import sys
sys.path.append(r"C:\Users\foivo\Documents\Coding\Scripts\fNIRS_fMRI_Project")
from utilities import fisher_z_transform, inverse_z_transform


def normalize_fmri_df(fmri_path, global_max):
    df = pd.read_excel(fmri_path, index_col=0)

    # fisher r to z the df
    df_z = fisher_z_transform(df)

    matrix = df_z.to_numpy()
    # remove the diagonal and fill the global max there instead
    np.fill_diagonal(matrix, global_max)

    # normalize the matrix
    normalized_matrix = matrix / global_max
    
    # create a new dataframe
    normalized_df = pd.DataFrame(normalized_matrix, index=df.index, columns=df.columns)
    return normalized_df

def find_max_fmri_value(fmri_path):
    df = pd.read_excel(fmri_path, index_col=0)

    # first we fisher z transform the dataframe to find the fisher z max because later we are going to normalize then in fisher z and then divide by the global fisher z max
    df_z = fisher_z_transform(df)
    # print(df_z.head())
    matrix = df_z.to_numpy()
    # remove the diagonal
    np.fill_diagonal(matrix, 0)

    # get the max absolute value
    max_value = np.max(np.abs(matrix))

    return max_value

def find_max_fnirs_value(fnirs_path):
    sheets = pd.ExcelFile(fnirs_path).sheet_names
    df_oxy = pd.read_excel(fnirs_path, index_col=0, sheet_name=sheets[0])

    # fisher r to z transform
    df_oxy_z = fisher_z_transform(df_oxy)

    matrix_oxy = df_oxy_z.to_numpy()
    # remove the diagonal
    np.fill_diagonal(matrix_oxy, 0)

    # get the max absolute value
    max_value_oxy = np.max(np.abs(matrix_oxy))

    df_deoxy = pd.read_excel(fnirs_path, index_col=0, sheet_name=sheets[1])

    # fisher r to z transform
    df_deoxy_z = fisher_z_transform(df_deoxy)

    matrix_deoxy = df_deoxy_z.to_numpy()

    # remove the diagonal
    np.fill_diagonal(matrix_deoxy, 0)
    max_value_deoxy = np.max(np.abs(matrix_deoxy))
    return max_value_oxy, max_value_deoxy

def normalize_fnirs_df(fnirs_path, global_max):
    sheets = pd.ExcelFile(fnirs_path).sheet_names
    df_oxy = pd.read_excel(fnirs_path, index_col=0, sheet_name=sheets[0])

    df_oxy_z = fisher_z_transform(df_oxy)

    matrix_oxy = df_oxy_z.to_numpy()
    # remove the diagonal and fill the global max there instead
    np.fill_diagonal(matrix_oxy, global_max)

    # normalize the matrix
    normalized_matrix_oxy = matrix_oxy / global_max
    
    # create a new dataframe
    normalized_oxy_df = pd.DataFrame(normalized_matrix_oxy, index=df_oxy.index, columns=df_oxy.columns)

    df_deoxy = pd.read_excel(fnirs_path, index_col=0, sheet_name=sheets[1])
    df_deoxy_z = fisher_z_transform(df_deoxy)

    matrix_deoxy = df_deoxy_z.to_numpy()
    # remove the diagonal and fill the global max there instead
    np.fill_diagonal(matrix_deoxy, global_max)

    # normalize the matrix
    normalized_matrix_deoxy = matrix_deoxy / global_max
    
    # create a new dataframe
    normalized_deoxy_df = pd.DataFrame(normalized_matrix_deoxy, index=df_deoxy.index, columns=df_deoxy.columns)

    return normalized_oxy_df, normalized_deoxy_df



# load the subject level folders of the matrices
fmri_master_folder_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\3_CorrelationMatrices"
fnirs_master_folder_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\1_CorrelationMatrices"

load_fmri_folders = [os.path.join(fmri_master_folder_path, f) for f in os.listdir(fmri_master_folder_path) if os.path.isdir(os.path.join(fmri_master_folder_path, f))]
load_fnirs_folders = [os.path.join(fnirs_master_folder_path, f) for f in os.listdir(fnirs_master_folder_path) if os.path.isdir(os.path.join(fnirs_master_folder_path, f))]
if len(load_fmri_folders) != 2 or len(load_fnirs_folders) !=2 or len(load_fmri_folders)!=len(load_fnirs_folders):
    raise ValueError("There should be exactly two folders for fMRI and two for fNIRS corresponding to Bivariate and Partial correlation matrices.")

fmri_par_global_max = 0
fmri_biv_global_max = 0
for fmri_folder in load_fmri_folders:
    fmri_files = [f for f in os.listdir(fmri_folder) if f.endswith('.xlsx')]
    print(f"Processing fMRI folder: {fmri_folder}")

    for fmri_file in fmri_files:
        fmri_matrix_path = os.path.join(fmri_folder, fmri_file)
        subject_max = find_max_fmri_value(fmri_matrix_path)

        if "Partial" in os.path.basename(fmri_folder):
            if subject_max >= fmri_par_global_max:
                fmri_par_global_max = subject_max
            # print(f"subject max: {subject_max}, current global max: {fmri_par_global_max}")
        if "Standard" in os.path.basename(fmri_folder):
            if subject_max >= fmri_biv_global_max:
                fmri_biv_global_max = subject_max
            # print(f"subject max: {subject_max}, current global max: {fmri_biv_global_max}")

print(f"fMRI Bivariate Global Max Abs Value: {fmri_biv_global_max}")
print(f"fMRI Partial Global Max Abs Value: {fmri_par_global_max}")

# now for fNIRS
fnirs_par_global_max_oxy = 0
fnirs_par_global_max_deoxy = 0
fnirs_biv_global_max_oxy = 0
fnirs_biv_global_max_deoxy = 0

for fnirs_folder in load_fnirs_folders:
    fnirs_files = [f for f in os.listdir(fnirs_folder) if f.endswith('.xlsx')]
    print(f"Processing fNIRS folder: {fnirs_folder}")
    for fnirs_file in fnirs_files:
        fnirs_matrix_path = os.path.join(fnirs_folder, fnirs_file)
        subject_max_oxy, subject_max_deoxy = find_max_fnirs_value(fnirs_matrix_path)
        if "Partial" in os.path.basename(fnirs_folder):
            if subject_max_oxy >= fnirs_par_global_max_oxy:
                fnirs_par_global_max_oxy = subject_max_oxy
            if subject_max_deoxy >= fnirs_par_global_max_deoxy:
                fnirs_par_global_max_deoxy = subject_max_deoxy
        if "Standard" in os.path.basename(fnirs_folder):
            if subject_max_oxy >= fnirs_biv_global_max_oxy:
                fnirs_biv_global_max_oxy = subject_max_oxy
            if subject_max_deoxy >= fnirs_biv_global_max_deoxy:
                fnirs_biv_global_max_deoxy = subject_max_deoxy

print(f"fNIRS Partial Global Max Abs Value HbO: {fnirs_par_global_max_oxy}")
print(f"fNIRS Partial Global Max Abs Value HbR: {fnirs_par_global_max_deoxy}")
print(f"fNIRS Bivariate Global Max Abs Value HbO: {fnirs_biv_global_max_oxy}")
print(f"fNIRS Bivariate Global Max Abs Value HbR: {fnirs_biv_global_max_deoxy}")
# now we have the global maximums for fMRI and fNIRS



fmri_out_master = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\3b_CorrelationMatrices_NormForGraph"
if not os.path.exists(fmri_out_master):
    os.makedirs(fmri_out_master)

for fmri_folder in load_fmri_folders:
    out_folder_name = os.path.join(fmri_out_master, os.path.basename(fmri_folder))
    if not os.path.exists(out_folder_name):
        os.makedirs(out_folder_name)
    
    fmri_files = [f for f in os.listdir(fmri_folder) if f.endswith('.xlsx')]
    for fmri_file in fmri_files:
        fmri_matrix_path = os.path.join(fmri_folder, fmri_file)
        normalized_df = pd.DataFrame()
        if "Partial" in os.path.basename(fmri_folder):
            normalized_df = normalize_fmri_df(fmri_matrix_path, fmri_par_global_max)

        if "Standard" in os.path.basename(fmri_folder):
            normalized_df = normalize_fmri_df(fmri_matrix_path, fmri_biv_global_max)
        out_path_file = os.path.join(out_folder_name, fmri_file.split('.')[0] + "_Normalized.xlsx")
        with pd.ExcelWriter(out_path_file) as writer:
            normalized_df.to_excel(writer, index=True)
        print(f"Saved normalized fMRI matrix to {out_path_file}")
        # quit()
fnirs_out_master = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\1b_CorrelationMatrices_NormForGraph"
if not os.path.exists(fnirs_out_master):
    os.makedirs(fnirs_out_master)

for fnirs_folder in load_fnirs_folders:
    out_folder_name = os.path.join(fnirs_out_master, os.path.basename(fnirs_folder))
    if not os.path.exists(out_folder_name):
        os.makedirs(out_folder_name)
    fnirs_files = [f for f in os.listdir(fnirs_folder) if f.endswith('.xlsx')]
    for fnirs_file in fnirs_files:
        fnirs_matrix_path = os.path.join(fnirs_folder, fnirs_file)
        normalized_oxy_df = pd.DataFrame()
        normalized_deoxy_df = pd.DataFrame()
        if "Partial" in os.path.basename(fnirs_folder):
            normalized_oxy_df, normalized_deoxy_df = normalize_fnirs_df(fnirs_matrix_path, fnirs_par_global_max_oxy)
        if "Standard" in os.path.basename(fnirs_folder):
            normalized_oxy_df, normalized_deoxy_df = normalize_fnirs_df(fnirs_matrix_path, fnirs_biv_global_max_oxy)
        out_path_file = os.path.join(out_folder_name, fnirs_file.split('.')[0] + "_Normalized.xlsx")
        with pd.ExcelWriter(out_path_file) as writer:
            normalized_oxy_df.to_excel(writer, index=True, sheet_name="HbO")
            normalized_deoxy_df.to_excel(writer, index=True, sheet_name="HbR")
print("Normalization Complete!")
print("New matrices saved")
