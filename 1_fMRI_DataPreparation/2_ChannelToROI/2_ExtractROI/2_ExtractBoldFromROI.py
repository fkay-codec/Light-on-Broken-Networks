"""
Script to extract BOLD time series from ROI spheres.
Loads combined ROI mask and fMRI data, extracts average BOLD signal 
from each ROI sphere, and saves results to Excel file.

NOTE to be scaled across multiple subjects
"""
from nilearn import plotting
from nilearn.image import load_img
from nilearn.maskers import NiftiLabelsMasker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


def extract_bold_save(roi_mask, functional_data):
    # Extract the subject name
    subj_id = os.path.basename(functional_data).split("_")[0] #! untested


    # Load the combined ROI file
    combined_img = load_img(roi_mask)

    # Load the functional data file
    func_img = load_img(functional_data)

    # Extract BOLD signals for each ROI
    print("     Extracting BOLD signals from ROIs...")
    masker = NiftiLabelsMasker(labels_img=combined_img, standardize=True)
    roi_time_series = masker.fit_transform(func_img)
    print("     Extraction complete.")
    # Get ROI labels
    roi_labels = masker.labels_

    # Create a DataFrame for easier handling
    roi_df = pd.DataFrame(roi_time_series, columns=[f'ROI_{label}' for label in roi_labels])

    # Create output directory if it doesn't exist
    output_folder = os.path.join("D:\Foivos\HCP_Data", "2_ExtractedBOLD")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # save the DataFrame as an excel file in which folder?
    excel_output = os.path.join(output_folder,f"{subj_id}_ROI_TimeSeries.xlsx")
    with pd.ExcelWriter(excel_output) as writer:
        roi_df.to_excel(writer, index=False)
    print("ROI time series saved to Excel file.")
    print(f"In folder: {output_folder}")

input_dir_roi_mask = r"D:\Foivos\HCP_Data\1_ROI_Sphere_Creation\Combined_ROIs.nii.gz"
input_dir_functional_data = r"D:\Foivos\HCP_Data\00_Subjects"

for file in os.listdir(input_dir_functional_data):
    if not file.endswith("_Rest3T.nii.gz"):
        continue  # Skip non-functional data files
    # call the function for each subject
    print(f"Processing file: {file}")
    functional_data = os.path.join(input_dir_functional_data, file)
    extract_bold_save(input_dir_roi_mask, functional_data)
    

