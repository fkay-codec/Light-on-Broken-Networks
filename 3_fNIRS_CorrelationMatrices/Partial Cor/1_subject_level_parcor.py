"""
Individual Subject Partial Correlation Matrix Pipeline for fNIRS Data

This script computes partial correlation matrices for individual subjects from preprocessed 
fNIRS SNIRF files. It processes both HbO and HbR signals separately to create subject-level 
connectivity matrices that will be used for downstream group-level statistical analysis.

Pipeline Overview:
1. Loads preprocessed SNIRF files from a selected folder
2. Extracts HbO and HbR channel data using MNE-Python
3. Removes short-distance channels (detector number > 28)
4. Computes partial correlation matrices using Pingouin for each hemoglobin type
5. Saves individual subject matrices to Excel files for group analysis

Partial Correlation:
For channels i and j, partial correlation PC(i,j) measures the correlation between 
channels i and j while controlling for the influence of all other channels in the dataset.
This removes confounding effects and reveals direct functional connectivity.

Output:
- Excel files named: {SubjectID}_ParCorM.xlsx
- Each file contains two sheets: 'HbO' and 'HbR'
- Matrices are symmetric with 1's on diagonal and partial correlation coefficients elsewhere

For theoretical background in preprocessing refer to MSc Thesis 'F.Kotsogiannis'.

NOTE: there was no thresholding based on p-values of the partial correlations, as this is not standard practice when we want to use graph metrics based on weighted matrices. However, this can be added if needed in future versions.
"""

import os
import time
import mne
import pingouin as pg
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import utils as ut


# Compute the partial correlation matrix
def compute_partial_correlation_matrix(df):
    """
    if i, j = 1, 2, ..., N are the N channels


    the partial correlation coef. between channels i and j irrespective of the effect N-2 
    PC(i, j) = correlation(i, j | all other channels) 

    To do that we need to calculate the covariance matrix (cov) of the data. 
    For N channels the PC will be equal to N(N-1)/2, which is obtained from the PC matrix   
    ----------
    here we use pingouin.partial_corr to calculate the partial correlation between each pair of channels
    
    Input:
    df: pandas dataframe with shape (n_samples, n_channels)
    Output:
    pc_df: pandas dataframe with shape (n_channels, n_channels) containing the partial correlation coefficients
        Only the pc_df is returned, because this is the one that we want. it is weighted
            
    """

    n_channels = len(df.columns)
    ch_names = df.columns.tolist()

    # We are going to build the pc_matrix: first the diagonal, then compute pc's, and then fill the symmetric side
    # save some computing time by calculating only half of the matrix (upper triangular), because if we are to perform this on many subjects across heme groups...
    pc_matrix = np.eye(n_channels)  # Initialize with identity matrix (1's on diagonal)

    for i in range(n_channels):
        for j in range(i + 1, n_channels):

            covariates = ch_names.copy()
            covariates.remove(ch_names[i])
            covariates.remove(ch_names[j])
            
            pc_result = pg.partial_corr(data=df, 
                                        x=ch_names[i], 
                                        y=ch_names[j], 
                                        covar=covariates)
            # print(f"Partial correlation between {ch_names[i]} and {ch_names[j]}: {pc_result['r'].values[0]}, p-value: {pc_result['p-val'].values[0]}")
            
            # here we can threshold it based on p-value if needed, but we dont currently.. still in testing phase.
            r_value = pc_result['r'].values[0]

            pc_matrix[i, j] = r_value
            pc_matrix[j, i] = r_value  # Symmetric entry


    # plot the matrices
    pc_df = pd.DataFrame(pc_matrix, index=ch_names, columns=ch_names)


    return pc_df

def extract_d_value(channel):
    """Designed to extract the d value from the channel name. Only works for channels with the format: D##_hbo or D##_hbr.
    keep as is and use a custom variation of this if needed. IMPORTANT DONT CHANGE THIS FUNCTION. JUST CREATE A NEW ONE IF NEEDED.
    """
    channel = channel.replace(" hbo", "").replace(" hbr", "")  # Remove the suffix to isolate the d value
    channel = channel.split('_')[-1]  # Extract the last part after the underscore, which is the: D##

    return channel[1:] # Extract the last part after the underscore, which is the d value

def detrend_polynomial(data):
    """Detrend data using polynomial fitting (2nd order)."""
    n_channels, n_samples = data.shape
    detrended_data = np.zeros_like(data, dtype=float)
    x = np.arange(n_samples)  # Time points
    for ch in range(n_channels):
        # Fit a 2nd order polynomial to the channel data
        coeffs = np.polyfit(x, data[ch], 2)
        # Evaluate the polynomial at the time points
        trend = np.polyval(coeffs, x)
        # Subtract the trend from the original data
        detrended_data[ch] = data[ch] - trend
    return detrended_data

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



def single_subject(input_file):
    """
    Compute the partial correlation matrix for a single subject fNIRS SNIRF file.
    
    Parameters:
    -----------
    input_file : str
        Path to the preprocessed SNIRF file for a single subject.
    
    Returns:
    --------
        Nothing. The function saves the partial correlation matrix to an Excel file.
        Output: Excel file with the partial correlation matrix for HbO and HbR in separate Sheets for single subject with name SubjectID_ParCorM.xlsx
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

    ## data could be detrended more if needed with the following line, not standard for par_cor, only for ICA; I believe it should be standard, but the differences are negligible; because it is already detrended for 1st polynominal during preprocessing
    #hbo_timeseries = detrend_polynomial(hbo_timeseries)  # Detrend the data (2nd order polynominal detrending)

    # create a dataframe with the timeseries and the channel names
    hbo_df = pd.DataFrame(hbo_timeseries.T, columns=hbo_channel_names)

    ## print what we did till now and check satori to verify
    # print(raw.times[-1])
    # print(hbo_df.head(), hbo_df.tail())
    # note everything looks good... continueing

    # Remove short channels before computing the partial correlation matrix
    hbo_df = remove_short_channels(hbo_df)

    # Keep only the channels that are in the fMRI analysis
    hbo_df = hbo_df[channels_in_fmri_list]

    hbo_partial_cor_matrix = compute_partial_correlation_matrix(hbo_df)

    # The same for HbR
    hbr_picks = mne.pick_types(raw.info, fnirs='hbr')
    hbr_channel_names = [raw.ch_names[pick] for pick in hbr_picks]
    hbr_timeseries = raw.get_data(picks=hbr_picks)

    hbr_df = pd.DataFrame(hbr_timeseries.T, columns=hbr_channel_names)
    hbr_df = remove_short_channels(hbr_df)
    hbr_df = hbr_df[channels_in_fmri_list]

    hbr_partial_cor_matrix = compute_partial_correlation_matrix(hbr_df)

    # Save the par. cor. matrix so we can manipulate it with other scripts if needed
    output_folder = r"D:\Foivos\fNIRS_fMRI_Study\fNIRS_data\1_CorrelationMatrices\2_PartialCor"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    with pd.ExcelWriter(os.path.join(output_folder, f"{subject_id}_PartialCorM.xlsx")) as writer:
        hbo_partial_cor_matrix.to_excel(writer, sheet_name='HbO')
        hbr_partial_cor_matrix.to_excel(writer, sheet_name='HbR')

# start time
# start= time.time()

# end = time.time()
# print("Partial Correlation Matrix computed and saved to Excel.")
# print(f"Time taken: {end - start} seconds")


input_folder =  r"C:\Users\foivo\Documents\Satori\SampleData\MULPA Dataset\3_Data_PrePro\01_RestingState_PrePro"

time_start = time.time()
for file in os.listdir(input_folder):
    if file.endswith(".snirf"):
        input_file = os.path.join(input_folder, file)
        # With this function we essential read the snirf file, extract the HbO and HbR channels, remove short channels, compute the par correlation, and then save it to an excel file.
        single_subject(input_file)
time_end = time.time()
print(f"Time taken for {os.path.basename(file)}: {time_end - time_start:.0f} seconds\n")