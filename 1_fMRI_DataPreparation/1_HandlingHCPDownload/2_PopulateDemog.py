"""When HCP data is downloaded/unzipped and reorganized with 1_HCP_FolderClean.py,
this script will read the HCP excel with the demo/behavioral data and populate another with only the age and gender of our downloaded subjects."""

import pandas as pd
import numpy as np
import os

# HCP excel currently in hard disc
HCP_excel = r"D:\Foivos\HCP_Data\Subjects\Subj_qualit_data.csv"

HCP_data = pd.read_csv(HCP_excel)
print("Available columns:", HCP_data.columns.tolist())

# fetch from the folder HCP_excel is located the list of the subject folders
HCP_folder = os.path.dirname(HCP_excel)
subj_folders = [f for f in os.listdir(HCP_folder) if not f.endswith(('.csv', '.xlsx', '.xls'))]
print("Subject folders found:", subj_folders)

subj_ids = [f.split('_')[0] for f in subj_folders]  # assuming folder names are like '102109_Rest3T'
print("Subject IDs:", subj_ids)

# Create a list to store demographic data (more efficient than append)
demog_list = []

for subj_id in subj_ids:
    subj_row = HCP_data[HCP_data['Subject'] == int(subj_id)]
    if not subj_row.empty:
        age = subj_row['Age'].values[0]
        gender = subj_row['Gender'].values[0]
        print(f"Subject {subj_id}: Age={age}, Gender={gender}")
        demog_list.append({'SubjectID': f"#{subj_id}", 'Age': age, 'Gender': gender})
    else:
        print(f"Warning: Subject {subj_id} not found in demographic data")

# Create DataFrame from list
demog_data = pd.DataFrame(demog_list)

# Save the demographic data to a new Excel file because I hate csvs
output_file = os.path.join(HCP_folder, 'selected_subjects_demographics.xlsx')
demog_data.to_excel(output_file, index=False)
print(f"\nDemographic data saved to: {output_file}")
print(f"Total subjects processed: {len(demog_data)}")