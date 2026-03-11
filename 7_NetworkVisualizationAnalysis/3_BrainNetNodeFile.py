"""
Here we parse the data to text files that BrainNet Viewer can read.
"""

import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import MinMaxScaler

def dataframe_to_brainnet(df, out, name):
    """
    Convert a DataFrame to BrainNet Viewer .node and .edge files.
    """
    df = df.copy()
    # create an output for the files to be saved
    out_folder = os.path.join(out, "2_BrainNetNodes")
    if not os.path.exists(out_folder):
        os.makedirs(out_folder, exist_ok=True)
    
    # First remap columns with the community assignment from 'A' 'B' 'C' to 1, 2, 3
    # we have in total 7 communities so we need 7 numbers
    community_map = {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7}
    gamma_list = []
    for col in df.columns:
        if "G=" not in col:
            continue
        df[col] = df[col].map(community_map)
        gamma_list.append(col)
    print(df.to_string())
    # drop any columns that are not needed for BrainNet
    df = df.drop('Pair Satori', axis=1, errors='ignore')
    # min max normalize the 'Nodal Strength Norm' column to be between 3 and 15 for better visualization
    # before that make the negatives zero
    df['Nodal Strength Norm'] = df['Nodal Strength Norm'].apply(lambda x: max(x, 0))

    scaler = MinMaxScaler(feature_range=(2.5, 5))
    df['Nodal Strength Norm'] = scaler.fit_transform(df[['Nodal Strength Norm']])

    # rename the cells in the column with the labels 
    df['ROI'] = [f.replace("_", "/") for f in df['ROI']]

    # print(df.head())
    
    # now we create the .node files, the number of .node depends on the number of columns with modules
    for gamma in gamma_list:
        print(gamma)
        clean_name = gamma.replace(".","p").replace("=","").replace(" ","").replace("|","_")
        node_file_name = name + "_" + clean_name + ".node"
        node_path = os.path.join(out_folder, node_file_name)
        print(f"Saving node file to: {node_path}")
        # lets create the .node file now
        with open(node_path, 'w') as f:
            for index, row in df.iterrows():
                # BrainNet Viewer .node format: x y z size color label
                x = row['X']
                y = row['Y']
                z = row['Z']
                size = round(row['Nodal Strength Norm'], 4)
                # size = 4.5
                color = row[gamma]
                label = row['ROI']
                # quit()
                f.write(f"{x}\t{y}\t{z}\t{color}\t{size}\t{label}\n")
        # quit()



fmri_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\7_GraphMetricsGroupLevel_Updated\BrainNetFiles\1_RawDataFrames"
fmri_master_dir = os.path.dirname(fmri_path)

fnirs_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\5_GraphMetricsGroupLevel_Updated\BrainNetFiles\1_RawDataFrames"
fnirs_master_dir = os.path.dirname(fnirs_path)

fmri_files = [f for f in os.listdir(fmri_path) if f.endswith('.xlsx')]
fnirs_files = [f for f in os.listdir(fnirs_path) if f.endswith('.xlsx')]

for file in fmri_files:
    if "Bivariate" in file:
        fmri_bivariate = pd.read_excel(os.path.join(fmri_path, file))
        dataframe_to_brainnet(fmri_bivariate, out=fmri_master_dir, name="fMRI_Bivariate")
    
    if "Partial" in file:
        fmri_partial = pd.read_excel(os.path.join(fmri_path, file))
        dataframe_to_brainnet(fmri_partial, out=fmri_master_dir, name="fMRI_Partial")


for file in fnirs_files:
    if "Bivariate" in file:
        fnirs_bivariate_hbo = pd.read_excel(os.path.join(fnirs_path, file), sheet_name="HbO_Bivariate_Corr")
        dataframe_to_brainnet(fnirs_bivariate_hbo, out=fnirs_master_dir, name="fNIRS_Bivariate_HbO")
        fnirs_bivariate_hbr = pd.read_excel(os.path.join(fnirs_path, file), sheet_name="HbR_Bivariate_Corr")
        dataframe_to_brainnet(fnirs_bivariate_hbr, out=fnirs_master_dir, name="fNIRS_Bivariate_HbR")
    
    if "Partial" in file:
        fnirs_partial_hbo = pd.read_excel(os.path.join(fnirs_path, file), sheet_name="HbO_Partial_Corr")
        dataframe_to_brainnet(fnirs_partial_hbo, out=fnirs_master_dir, name="fNIRS_Partial_HbO")
        fnirs_partial_hbr = pd.read_excel(os.path.join(fnirs_path, file), sheet_name="HbR_Partial_Corr")
        dataframe_to_brainnet(fnirs_partial_hbr, out=fnirs_master_dir, name="fNIRS_Partial_HbR")
