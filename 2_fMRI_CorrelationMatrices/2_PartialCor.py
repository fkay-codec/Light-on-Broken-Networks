""""
In this script a partial correlation matrix will be computed
Input: 
    Folder containing the timeseries fMRI data
    Pandas dataframe with each column as ROI timeseries
Output:
    Excel file with the correlation matrix ROI x ROI
"""

from rich.progress import Progress
import time
import os
import pandas as pd
import numpy as np
import pingouin as pg

def compute_partial_correlation_matrix(df):
    """
    Computes the partial correlation matrix for all variable pairs in a dataframe.
    Columns are variables, rows are datapoints.
    
    Parameters:
        df (pd.DataFrame): DataFrame with shape (n_samples, n_variables)
    
    Returns:
        pc_df (pd.DataFrame): Partial correlation matrix (n_variables x n_variables)
    """
    n_vars = len(df.columns)
    var_names = df.columns.tolist()
    pc_matrix = np.eye(n_vars)  # Diagonal = 1

    for i in range(n_vars):
        for j in range(i + 1, n_vars):
            covariates = var_names.copy()
            covariates.remove(var_names[i])
            covariates.remove(var_names[j])
            pc_result = pg.partial_corr(data=df, x=var_names[i], y=var_names[j], covar=covariates)
            r_value = pc_result['r'].values[0]
            pc_matrix[i, j] = r_value
            pc_matrix[j, i] = r_value  # Symmetric

    pc_df = pd.DataFrame(pc_matrix, index=var_names, columns=var_names)
    return pc_df

def clean_fmri_dataframe(df):
    """
    needs dataframe to work [rows: time points, columns: ROIs]
    1. convert names with ROI1.0 to ROI1
    2. add corresponding fnirs channel name on the column names
    """
    # clean the column names
    df.columns = [col.replace(".0", "").replace("_", "") for col in df.columns]
    # add the corresponding fnirs channel name by itterating through the folder that contains the fMRI masks created from fNRIRS channels. In this folder the file names contain roi-channel names
    dir = r"D:\Foivos\fNIRS_fMRI_Study\fMRI_data\1_ROI_Sphere_Creation"
    for col in df.columns:
        for file in os.listdir(dir):
            if not file.endswith("_point.nii.gz"):
                continue  # Skip non-functional data files
            file_roi = file.split("_")[0]  # Extract the ROI part from the filename
            if col == file_roi:
                # print(f"Matching {col} with {file}")
                # rename the whole column now with the filename - the irrelavant parts
                new_col_name = file.replace("_point.nii.gz", "").replace("-","")
                # print(f"Renaming {col} to {new_col_name}")
                df = df.rename(columns={col: new_col_name})
    return df 


def single_subject_analysis(input_file):
    """
    Compute the ROI wise partial correlation matrix for a single subject. 
        Maybe add the equivelant channel beside the ROI name for easier identification/comparability?
    
    Parameters:
    -----------
    input_file : str
        Path to the preprocessed SNIRF file for a single subject.
    
    Returns:
    --------
        Nothing. The function saves the correlation matrix to an Excel file.
    """

    subject_id = os.path.basename(input_file).split('_')[0]

    # load the dataframe
    df = pd.read_excel(input_file)
    # print(df.head())

    # with this function, i clean the fmri dataframe by ranaming the columns to match the fnirs channel names. eg., ROI1_S1D1               
    df = clean_fmri_dataframe(df)
    # compute the roi wise correlation
    df_corr = compute_partial_correlation_matrix(df)

    # Save the cor. matrix so we can manipulate it with other scripts if needed
    output_folder = r"D:\Foivos\fNIRS_fMRI_Study\fMRI_data\3_CorrelationMatrices"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_ofcor = os.path.join(output_folder, "2_PartialCor")
    if not os.path.exists(output_ofcor):
        os.makedirs(output_ofcor)
    with pd.ExcelWriter(os.path.join(output_ofcor, f"{subject_id}_PartialCorM.xlsx")) as writer:   
        df_corr.to_excel(writer, sheet_name='PartialCorM')

    print(f"    Correlation matrix saved for subject {subject_id} in folder: {output_ofcor}")


# folder that contains the excel files with the timeseries data
input_folder = r"D:\Foivos\fNIRS_fMRI_Study\fMRI_data\2_ExtractedBOLD"
files = [f for f in os.listdir(input_folder) if f.endswith(".xlsx")]
with Progress() as progress:
    task  = progress.add_task(f"Processing ...", total=len(files))
    for file in os.listdir(input_folder):
        if file.endswith(".xlsx"):
            input_file = os.path.join(input_folder, file)
            print(f"Processing file: {input_file}")
            single_subject_analysis(input_file)
            progress.update(task, advance=1)