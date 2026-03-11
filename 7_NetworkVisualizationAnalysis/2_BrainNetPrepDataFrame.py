"""
The goal of this script is to create a template for BrainNet Viewer

We want to vizualize the graph in terms of modules and nodal strength

To do that we need a compatible .node file for BrainNet to read and create

Thus we need:
1. Read the excel with the channel names and coordinates
2. Add the ROI numbers in the dataframe that match the channel names based on our scheme
3. Read the nodal strength values and add them to the dataframe
4. Read the modularity values and add them to the dataframe
5. Save the dataframe as a .node file

*In this specific script we are creating a master dataframe with all relevant information for the creation of the .node
*Channel names | X | Y | Z | ROI number/Channel | Modularity Column 1 | Modularity Column 2 | ... | Nodal Strength Normalized
Note: the reasoning stems from the fact that when we visualize both fnirs and fmri the Nodal Strength Normalized will be the same if we scale it to 1. So for later .node files we need all the data to normalize fNIRS/fMRI node str with the same reference
"""

import os
import numpy as np
import pandas as pd
import sys

# lets test first
def prepare_dataframe():

    """
    Prepare the master dataframe.
    the output will be a dataframe with columns: Pair Satori, X, Y, Z, ROI number/Channel
    """
    # this is the excel file that has the channel names and coordinates in MNI space, extracted from the fOLD toolbox and then dropped the rows(channels) that dont contain any MNI coords.
    channel_to_brain_df = pd.read_excel(r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\Channels to ROI\Channels to Brain Areas using fOLD_updated_droppedmissing.xlsx")

    # this contains the folder directory in which we assigned channels a specific ROI number
    roi_dir = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\1_ROI_Sphere_Creation"
    roi_files = [f for f in os.listdir(roi_dir) if "_point" in f]
    cleaned_name = [f.replace("_point.nii.gz", "") for f in roi_files]
    rois = [f.split("_")[0] for f in cleaned_name]
    channels = [f.split("_")[1] for f in cleaned_name]

    # hard check
    for clean, roi, ch in zip(cleaned_name, rois, channels):
        if clean != f"{roi}_{ch}":
            print(f"Filename {clean} does not match ROI {roi} and Channel {ch}")
            print("         Catastrophic failure")
            quit()
    # successfull mapping
    channels_to_roi = dict(zip(rois, channels))
    # hard check of dictionary
    for key, ch_orig in zip(channels_to_roi.items(), cleaned_name):
        if f"{key[0]}_{key[1]}" != ch_orig:
            print(f"Dictionary key-value pair {key[0]}_{key[1]} does not match filename {ch_orig}")
            print("         Catastrophic failure")
            quit()
    # quit

    print("Successfully mapped ROI files to channels")


    master_df = channel_to_brain_df[["Pair Satori", "X", "Y", "Z"]].copy()
    print(master_df.head())
    # now if you find key in Pair_Satori, assign ROI there
    master_df["ROI"] = pd.Series(dtype=str)
    for sat_pair in master_df["Pair Satori"]:
        if sat_pair not in channels_to_roi.values():
            print(f"Channel {sat_pair} not found in ROI mapping")
            print("         Catastrophic failure")
            quit()
        for roi, channel in channels_to_roi.items():
            # print(roi, channel)
            if channel == sat_pair:
                # print("found it", channel, roi, sat_pair)
                row_index = master_df[master_df["Pair Satori"] == sat_pair].index[0]
                master_df.loc[row_index,"ROI"] = f"{roi}_{channel}"
    # print(master_df.head())
    #! Hard Check
    #* The Hard Checks although are repetitive, i want to be 100000% sure that everything is mapped correctly, it is really crucial
    for sat_col, roi_col in zip(master_df["Pair Satori"], master_df["ROI"]):
        roi_col_clean = roi_col.split("_")[1]  # get the channel part
        if sat_col != roi_col_clean:
            print(f"Channel {sat_col} does not match ROI column {roi_col}")
            print("         Catastrophic failure")
            quit()
    print("Successfully assigned ROI numbers to channels")

    master_df["ROI"] = master_df["ROI"].str.replace("-", "") 
    return master_df

def check_compatibility_of_list(list1, list2):
    """
    Check if the channel names in the list1 are compatible with the ones in the list2, in terms of order
    """
    for chan1, chan2 in zip(list1, list2):

        if chan1 != chan2:
            print(f"Channel {chan1} is not compatible with {chan2}")
            print("         Catastrophic failure")
            quit()
    print("Lists are compatible")
    return True

def check_compatibility(df1, df2):
    """
    Check if the channel names in the df are compatible with the ones in the modularity df, in terms of order
    """
    df1 = df1.copy()
    df2 = df2.copy()
    df1["Pair Satori"] = df1["Pair Satori"].str.replace("-", "").str.replace("_", "")  # remove any - or _ from the channel names

    if "ROI" in df2.index[0]:
        df2.index = df2.index.str.split("_").str[1]  # remove the ROI part

    for chan1, chan2 in zip(df1["Pair Satori"], df2.index):

        if chan1 != chan2:
            print(f"Channel {chan1} is not compatible with {chan2}")
            print("         Catastrophic failure")
            quit()
    print("Dataframes are compatible")
    return True

# Add modularity values in the dataframe

# to do that we have to open the fMRI excel file, one for partial, one for bivariate.
# then the second sheet has the modularity per gamma...
# we are selecting gamma 0.4-0.6 for partial || 0
def add_modularity_to_df(df, modul_df):
    """
    Check how many unique communities we have in the column... we want the ones that have between 5-7

    If we find such a column, we add it to the dataframe and return the new dataframe
    """
    check_compatibility(df, modul_df)
    df = df.copy()

    for col in modul_df.columns:
        col_to_add = None
        unique_com = modul_df[col].nunique()
        if 5 <= unique_com <= 7:
            # print(f"Column {col} has {unique_com} unique communities")

            # add the selected col to the dataframe
            col_to_add = modul_df[col].reset_index(drop=True)
            df = pd.concat([df, col_to_add], axis=1)
                                      
    return df

def process_nodal_strength(df):
    """
    Process the dataframe

    Extract Nodal Str and then normalize it.

    At the end return a series with the normalized values
    
    """
    graph_df = df.copy()

    nodal_str = graph_df["Normalized Nodal Strength"]

    return nodal_str

# now that we have everything in the dataframe we can create the .node files with the relevant information
def create_node_file(df, out_save):
    """
    Create BrainNet compatible .node files based on the dataframe provided

    Number of files created = number of modularity columns in the dataframe
    out_save = directory to save the files

    Save in out_save folder

    .node file format:
    x y z color size label
    NOTE: instead of 'space' use 'tab' to separate the columns
    
    Consistent coloring scheme:
    'A' = 1; 'B' = 2; 'C' = 3; 'D' = 4; 'E' = 5; 'F' = 6; 'G' = 7
    """



df = prepare_dataframe()

print("Processing fMRI files")
fmri_path_modularity =r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\7_GraphMetricsGroupLevel_Updated"
fmri_modul_files=[f for f in os.listdir(fmri_path_modularity) if f.endswith(".xlsx")]
print(fmri_modul_files)
for file in fmri_modul_files:

    if "partial" in file.lower():
        print("     Processing: Partial Correlation Modularity")
        fmri_par_path = os.path.join(fmri_path_modularity, file)
        fmri_par_mod_df = pd.read_excel(fmri_par_path, sheet_name=1, index_col=0)
        fmri_par_template_df = df.copy()

        fmri_par_df = add_modularity_to_df(fmri_par_template_df, fmri_par_mod_df)
        fmri_par_df_roi_list = fmri_par_df["ROI"].tolist()


        fmri_par_graph_df = pd.read_excel(fmri_par_path, sheet_name=0, index_col=0)
        fmri_par_graph_df_index_list = fmri_par_graph_df.index.tolist()
        check_compatibility_of_list(fmri_par_df_roi_list, fmri_par_graph_df_index_list)


        fmri_par_nodal_df = process_nodal_strength(fmri_par_graph_df)
        fmri_par_nodal_df_index_list = fmri_par_nodal_df.index.tolist()
        check_compatibility_of_list(fmri_par_df_roi_list, fmri_par_nodal_df_index_list)

        fmri_par_nodal_df.reset_index(drop=True, inplace=True)

        fmri_par_df["Nodal Strength Norm"] = fmri_par_nodal_df

    if "standard" in file.lower():
        fmri_biv_path = os.path.join(fmri_path_modularity, file)
        fmri_biv_mod_df = pd.read_excel(fmri_biv_path, sheet_name=1, index_col=0)
        print("Processing: Bivariate Correlation Modularity")
        fmri_biv_template_df = df.copy()
        
        fmri_biv_df = add_modularity_to_df(fmri_biv_template_df, fmri_biv_mod_df)
        fmri_biv_df_roi_list = fmri_biv_df["ROI"].tolist()

        fmri_biv_graph_df = pd.read_excel(fmri_biv_path, sheet_name=0, index_col=0)
        fmri_biv_graph_df_index_list = fmri_biv_graph_df.index.tolist()
        check_compatibility_of_list(fmri_biv_df_roi_list, fmri_biv_graph_df_index_list)

        fmri_biv_nodal_df = process_nodal_strength(fmri_biv_graph_df)
        fmri_biv_nodal_df_index_list = fmri_biv_nodal_df.index.tolist()
        check_compatibility_of_list(fmri_biv_df_roi_list, fmri_biv_nodal_df_index_list)
        fmri_biv_nodal_df.reset_index(drop=True, inplace=True)

        fmri_biv_df["Nodal Strength Norm"] = fmri_biv_nodal_df

print("Finished processing fMRI files\n\n")
# print(fmri_par_df.head())
# print(fmri_biv_df.head())


# now the same for fNRIS
fnirs_modul_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\5_GraphMetricsGroupLevel_Updated"
fnirs_modul_files = [f for f in os.listdir(fnirs_modul_path) if f.endswith(".xlsx")]
print("Processing fNIRS files:")
for file in fnirs_modul_files:
    if "partial" in file.lower():
        print("     Processing: Partial Correlation Modularity")
        # HbO
        fnirs_par_path = os.path.join(fnirs_modul_path, file)
        fnirs_par_hbo_mod_df = pd.read_excel(fnirs_par_path, sheet_name="HbO_Modularity", index_col=0)
        fnirs_par_hbo_mod_df.index = fnirs_par_hbo_mod_df.index.str.replace("_", "") #.str.replace("_", "")  # remove any - or _ from the channel names
        fnirs_par_hbo_template_df = df.copy()

        fnirs_par_hbo_df = add_modularity_to_df(fnirs_par_hbo_template_df, fnirs_par_hbo_mod_df)
        fnirs_par_hbo_df_roi_list = fnirs_par_hbo_df["ROI"].tolist()
        fnirs_par_hbo_df_roi_list_comp = [f.split("_")[1] for f in fnirs_par_hbo_df_roi_list]  # remove ROI part if exists

        fnirs_par_hbo_graph_df = pd.read_excel(fnirs_par_path, sheet_name="HbO_GraphMetrics", index_col=0)
        fnirs_par_hbo_graph_df.index = fnirs_par_hbo_graph_df.index.str.replace("_", "") #.str.replace("_", "")  # remove any - or _ from the channel names
        fnirs_par_hbo_graph_df_index_list = fnirs_par_hbo_graph_df.index.tolist()
        check_compatibility_of_list(fnirs_par_hbo_df_roi_list_comp, fnirs_par_hbo_graph_df_index_list)

        fnirs_par_hbo_nodal_df = process_nodal_strength(fnirs_par_hbo_graph_df)
        fnirs_par_hbo_nodal_df_index_list = fnirs_par_hbo_nodal_df.index.tolist()
        check_compatibility_of_list(fnirs_par_hbo_df_roi_list_comp, fnirs_par_hbo_nodal_df_index_list)
        # print(fnirs_par_hbo_nodal_df.head()
        fnirs_par_hbo_nodal_df.reset_index(drop=True, inplace=True)
        fnirs_par_hbo_df["Nodal Strength Norm"] = fnirs_par_hbo_nodal_df
    
        # HbR
        fnirs_par_hbr_mod_df = pd.read_excel(fnirs_par_path, sheet_name="HbR_Modularity", index_col=0)
        fnirs_par_hbr_mod_df.index = fnirs_par_hbr_mod_df.index.str.replace("_", "") #.str.replace("_", "")  # remove any - or _ from the channel names
        fnirs_par_hbr_template_df = df.copy()   
        fnirs_par_hbr_df = add_modularity_to_df(fnirs_par_hbr_template_df, fnirs_par_hbr_mod_df)
        fnirs_par_hbr_df_roi_list = fnirs_par_hbr_df["ROI"].tolist()
        fnirs_par_hbr_df_roi_list_comp = [f.split("_")[1] for f in fnirs_par_hbr_df_roi_list]  # remove ROI part if exists
        fnirs_par_hbr_graph_df = pd.read_excel(fnirs_par_path, sheet_name="HbR_GraphMetrics", index_col=0)
        fnirs_par_hbr_graph_df.index = fnirs_par_hbr_graph_df.index.str.replace("_", "") #.str.replace("_", "")  # remove any - or _ from the channel names
        fnirs_par_hbr_graph_df_index_list = fnirs_par_hbr_graph_df.index.tolist()
        check_compatibility_of_list(fnirs_par_hbr_df_roi_list_comp, fnirs_par_hbr_graph_df_index_list)
        fnirs_par_hbr_nodal_df = process_nodal_strength(fnirs_par_hbr_graph_df)
        fnirs_par_hbr_nodal_df_index_list = fnirs_par_hbr_nodal_df.index.tolist()
        check_compatibility_of_list(fnirs_par_hbr_df_roi_list_comp, fnirs_par_hbr_nodal_df_index_list)
        fnirs_par_hbr_nodal_df.reset_index(drop=True, inplace=True) 
        fnirs_par_hbr_df["Nodal Strength Norm"] = fnirs_par_hbr_nodal_df    
        
        print(fnirs_par_hbo_df.head())
        print(fnirs_par_hbr_df.head())
    if "standard" in file.lower():
        print("     Processing: Bivariate Correlation Modularity")
        fnirs_biv_path = os.path.join(fnirs_modul_path, file)
        # HbO
        fnirs_biv_hbo_mod_df = pd.read_excel(fnirs_biv_path, sheet_name="HbO_Modularity", index_col=0)
        fnirs_biv_hbo_mod_df.index = fnirs_biv_hbo_mod_df.index.str.replace("_", "") #.str.replace("_", "")  # remove any - or _ from the channel names
        fnirs_biv_hbo_template_df = df.copy()
        fnirs_biv_hbo_df = add_modularity_to_df(fnirs_biv_hbo_template_df, fnirs_biv_hbo_mod_df)
        fnirs_biv_hbo_df_roi_list = fnirs_biv_hbo_df["ROI"].tolist()
        fnirs_biv_hbo_df_roi_list_comp = [f.split("_")[1]   for f in fnirs_biv_hbo_df_roi_list]  # remove ROI part if exists
        fnirs_biv_hbo_graph_df = pd.read_excel(fnirs_biv_path, sheet_name="HbO_GraphMetrics", index_col=0)
        fnirs_biv_hbo_graph_df.index = fnirs_biv_hbo_graph_df.index.str.replace("_", "") #.str.replace("_", "")  # remove any - or _ from the channel names
        fnirs_biv_hbo_graph_df_index_list = fnirs_biv_hbo_graph_df.index.tolist()
        check_compatibility_of_list(fnirs_biv_hbo_df_roi_list_comp, fnirs_biv_hbo_graph_df_index_list)
        fnirs_biv_hbo_nodal_df = process_nodal_strength(fnirs_biv_hbo_graph_df)
        fnirs_biv_hbo_nodal_df_index_list = fnirs_biv_hbo_nodal_df.index.tolist()
        check_compatibility_of_list(fnirs_biv_hbo_df_roi_list_comp, fnirs_biv_hbo_nodal_df_index_list)
        fnirs_biv_hbo_nodal_df.reset_index(drop=True, inplace=True)    
        fnirs_biv_hbo_df["Nodal Strength Norm"] = fnirs_biv_hbo_nodal_df
        # HbR   
        fnirs_biv_hbr_mod_df = pd.read_excel(fnirs_biv_path, sheet_name="HbR_Modularity", index_col=0)
        fnirs_biv_hbr_mod_df.index = fnirs_biv_hbr_mod_df.index.str.replace("_", "") #.str.replace("_", "")  # remove any - or _ from the channel names
        fnirs_biv_hbr_template_df = df.copy()
        fnirs_biv_hbr_df = add_modularity_to_df(fnirs_biv_hbr_template_df, fnirs_biv_hbr_mod_df)
        fnirs_biv_hbr_df_roi_list = fnirs_biv_hbr_df["ROI"].tolist()
        fnirs_biv_hbr_df_roi_list_comp = [f.split("_")[1]   for f in fnirs_biv_hbr_df_roi_list]  # remove ROI part if exists
        fnirs_biv_hbr_graph_df = pd.read_excel(fnirs_biv_path, sheet_name="HbR_GraphMetrics", index_col=0)
        fnirs_biv_hbr_graph_df.index = fnirs_biv_hbr_graph_df.index.str.replace("_", "") #.str.replace("_", "")  # remove any - or _ from the channel names
        fnirs_biv_hbr_graph_df_index_list = fnirs_biv_hbr_graph_df.index.tolist()
        check_compatibility_of_list(fnirs_biv_hbr_df_roi_list_comp, fnirs_biv_hbr_graph_df_index_list)
        fnirs_biv_hbr_nodal_df = process_nodal_strength(fnirs_biv_hbr_graph_df)
        fnirs_biv_hbr_nodal_df_index_list = fnirs_biv_hbr_nodal_df.index.tolist()
        check_compatibility_of_list(fnirs_biv_hbr_df_roi_list_comp, fnirs_biv_hbr_nodal_df_index_list)
        fnirs_biv_hbr_nodal_df.reset_index(drop=True, inplace=True)
        fnirs_biv_hbr_df["Nodal Strength Norm"] = fnirs_biv_hbr_nodal_df

print("Finished processing fNIRS files\n\n")

print(fmri_par_df.head())
print(fmri_biv_df.head())
print(fnirs_par_hbo_df.head())
print(fnirs_par_hbr_df.head())
print(fnirs_biv_hbo_df.head())
print(fnirs_biv_hbr_df.head())

# save everything in excel

#* fMRI save directory
out_save_dir = os.path.join(fmri_path_modularity, "BrainNetFiles", "1_RawDataFrames")
if not os.path.exists(out_save_dir):
    os.makedirs(out_save_dir)
# save the dataframes as excel files for future reference
fmri_par_df.to_excel(os.path.join(out_save_dir, "fMRI_Partial_Corr_Master_DF.xlsx"), index=False)
fmri_biv_df.to_excel(os.path.join(out_save_dir, "fMRI_Bivariate_Corr_Master_DF.xlsx"), index=False)


#* fNIRS save the dataframes as excel files for future reference
# Save directory
out_save_dir = os.path.join(fnirs_modul_path, "BrainNetFiles", "1_RawDataFrames")
if not os.path.exists(out_save_dir):
    os.makedirs(out_save_dir)
with pd.ExcelWriter(os.path.join(out_save_dir, "fNIRS_Partial_Corr_Master_DF.xlsx")) as writer:
    fnirs_par_hbo_df.to_excel(writer, sheet_name="HbO_Partial_Corr", index=False)
    fnirs_par_hbr_df.to_excel(writer, sheet_name="HbR_Partial_Corr", index=False)
with pd.ExcelWriter(os.path.join(out_save_dir, "fNIRS_Bivariate_Corr_Master_DF.xlsx")) as writer:
    fnirs_biv_hbo_df.to_excel(writer, sheet_name="HbO_Bivariate_Corr", index=False)
    fnirs_biv_hbr_df.to_excel(writer, sheet_name="HbR_Bivariate_Corr", index=False)

print("Dataframes saved as excel files")

