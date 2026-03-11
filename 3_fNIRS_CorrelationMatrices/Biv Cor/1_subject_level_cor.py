""""This script is similar to the SBA_Corr #1a CrossCorrelationPerSubject.py
#! RESULTS ARE VERIFIED WITH PREVIOUS CROSS CORRELATION RESULTS FROM MSc THESIS F.K.
see SBA-COR-SCRIPT: 1a CrossCorrelationPerSubject.py
"""


import sys
from PySide6.QtWidgets import QApplication, QFileDialog
import os
import time
import mne
import pandas as pd
import numpy as np


def select_folder_path(prompt):
    # Open dialog to select the folder of interest and return its path
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)

    # Open a dialog to select a folder 
    folder_path = QFileDialog.getExistingDirectory(None, prompt, r"C:\Users\foivo\Documents\Satori\SampleData\Tilburg_Paper")
    return folder_path

def extract_d_value(channel):
    """Designed to extract the d value from the channel name. Only works for channels with the format: D##_hbo or D##_hbr.
    keep as is and use a custom variation of this if needed. IMPORTANT DONT CHANGE THIS FUNCTION. JUST CREATE A NEW ONE IF NEEDED.
    """
    channel = channel.replace(" hbo", "").replace(" hbr", "")  # Remove the suffix to isolate the d value
    channel = channel.split('_')[-1]  # Extract the last part after the underscore, which is the: D##

    return channel[1:] # Extract the last part after the underscore, which is the d value

# #def detrend_polynomial(data):
#     """Detrend data using polynomial fitting (2nd order)."""
#     n_channels, n_samples = data.shape
#     detrended_data = np.zeros_like(data, dtype=float)
#     x = np.arange(n_samples)  # Time points
#     for ch in range(n_channels):
#         # Fit a 2nd order polynomial to the channel data
#         coeffs = np.polyfit(x, data[ch], 2)
#         # Evaluate the polynomial at the time points
#         trend = np.polyval(coeffs, x)
#         # Subtract the trend from the original data
#         detrended_data[ch] = data[ch] - trend
#     return detrended_data

def remove_short_channels(df):
    """
    needs dataframe to work [rows: time points, columns: channels]

    1. remove 'hbo' or 'hbr' from channel names
    2. identify short channels
    3. drop short channels from dataframe
    4. return the cleaned dataframe    
    """
    # remove hbo from channel names
    df.columns = [col.replace(" hbo", "").replace(" hbr", "") for col in df.columns]

    columns_to_drop = []
    # identify short channels (d > 28)
    for col in df.columns:
        channel = col.split('_')[-1]  # Extract the last part after the underscore, which is the:
        # print(channel)
        d_value = int(channel[1:])  # Extract the d value
        if d_value > 28:
            columns_to_drop.append(col)
    # drop short channels from dataframe
    df = df.drop(columns=columns_to_drop)
    return df

def single_subject_cross_corr(input_file):
    """
    Compute the channel wise correlation matrix for a single subject fNIRS SNIRF file.
    
    Parameters:
    -----------
    input_file : str
        Path to the preprocessed SNIRF file for a single subject.
    
    Returns:
    --------
        Nothing. The function saves the correlation matrix to an Excel file.
        Output: Excel file with the correlation matrix for HbO and HbR in separate Sheets for single subject with name SubjectID_CorM.xlsx
    """
    #! To make the results comparable with the fMRI we have to drop the channels that are not in the fMRI analysis
    channels_in_fmri=r"D:\Foivos\fNIRS_fMRI_Study\Channels to Brain Areas using fOLD_updated_droppedmissing.xlsx"
    channels_in_fmri_df = pd.read_excel(channels_in_fmri)
    channels_in_fmri_list = channels_in_fmri_df['Pair Satori'].tolist()
    # replace "-" with "" in the list
    channels_in_fmri_list = [ch.replace("-", "_") for ch in channels_in_fmri_list]



    print(f"Processing file: {input_file}")
    # I am extracting the subj ID in order to write the new file and have an indentifiable naming convention
    subject_id = os.path.basename(input_file).split('_')[0]
    print(f"Subject ID: {subject_id}")

    # read the .snirf file with MNE
    raw = mne.io.read_raw_snirf(input_file, preload=True)

    # extract the HbO signals from the channels we want and their names
    hbo_picks = mne.pick_types(raw.info, fnirs='hbo')
    hbo_channel_names = [raw.ch_names[pick] for pick in hbo_picks]
    hbo_timeseries = raw.get_data(picks=hbo_picks)

    ## data could be detrended more if needed with the following line, not standard for cor, only for ICA; I believe it should be standard, but the differences are negligible; because it is already detrended for 1st polynominal during preprocessing
    #hbo_timeseries = detrend_polynomial(hbo_timeseries)  # Detrend the data (2nd order polynominal detrending)

    # create a dataframe with the timeseries and the channel names
    hbo_df = pd.DataFrame(hbo_timeseries.T, columns=hbo_channel_names)

    ## print what we did till now and check satori to verify

    # note everything looks good... continueing

    # Remove short channels before computing the  correlation matrix
    hbo_df = remove_short_channels(hbo_df)

    # Keep only the channels that are in the fMRI analysis
    hbo_df = hbo_df[channels_in_fmri_list]


    # compute the channel wise correlation
    hbo_corr_df = hbo_df.corr()

    # The same for HbR
    hbr_picks = mne.pick_types(raw.info, fnirs='hbr')
    hbr_channel_names = [raw.ch_names[pick] for pick in hbr_picks]
    hbr_timeseries = raw.get_data(picks=hbr_picks)

    hbr_df = pd.DataFrame(hbr_timeseries.T, columns=hbr_channel_names)
    hbr_df = remove_short_channels(hbr_df)
    hbr_df = hbr_df[channels_in_fmri_list]
    hbr_corr_df = hbr_df.corr()

    # Save the cor. matrix so we can manipulate it with other scripts if needed
    output_folder = r"D:\Foivos\fNIRS_fMRI_Study\fNIRS_data\1_CorrelationMatrices\1_StandardCor"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    with pd.ExcelWriter(os.path.join(output_folder, f"{subject_id}_StandardCorM.xlsx")) as writer:
        hbo_corr_df.to_excel(writer, sheet_name='HbO')
        hbr_corr_df.to_excel(writer, sheet_name='HbR')


input_folder = r"C:\Users\foivo\Documents\Satori\SampleData\MULPA Dataset\3_Data_PrePro\01_RestingState_PrePro"
# select_folder_path("Select the folder containing preprocessed SNIRF files")
time_start = time.time()

for file in os.listdir(input_folder):
    if file.endswith(".snirf"):
        input_file = os.path.join(input_folder, file)
        # With this function we essential read the snirf file, extract the HbO and HbR channels, remove short channels, compute the correlation, and then save it to an excel file.
        single_subject_cross_corr(input_file)
time_end = time.time()
print(f"Time taken for {os.path.basename(file)}: {time_end - time_start:.0f} seconds\n")