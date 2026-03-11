"""In this scrip we will extract from the HCP folders the relevant files and copy them to a new folder structure for easier handling
HCP folder structure with RS fMRI 3T Data:

SubjID_Rest3TRecommended
    /SubjID
        /MNINonLinear
            /fsaverage_LR32k // contains the surface files/folders // #! NOT NEEDED
            /Results
                /rfMRI_REST
                /rfMRI_REST1_LR // #* target folder for LP run from session 1
                    ...
                    /rfMRI_REST1_LR_hp2000_clean_rclean_tclean.nii.gz #* target NIfTI file for LP run from session 1
                    ...
                /rfMRI_REST1_RL
                /rfMRI_REST2_LR
                /rfMRI_REST2_RL

                
Output folder Structure:
Subjects
    SubjID_NIfTIData.nii.gz #* with new naming convention for easier handling. EG: SubjID_Rest3T.nii.gz
    ...


#! When the script completes, verify the output, delete the input and transfer the files to a storage disc for safekeeping
"""

import os
import shutil

# Define source and target directories

input_dir = r"D:\Foivos\HCP_Data\RAW"
# r"C:\Users\foivo\Desktop\Ubuntu Sharing\fMRI data\Subjects_RawHCPDownload" # this is where the original HCP folders are stored

output_dir = r"D:\Foivos\HCP_Data\Subjects"
# r"C:\Users\foivo\Desktop\Ubuntu Sharing\fMRI data\Subjects"  # this is where we are going to transfer our files

# List all subject directories in the input directory
subject_dirs = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
print(f"Found {len(subject_dirs)} subject directories.")

# Iterate through each subject directory
for folder in subject_dirs:
    if folder.endswith(('.zip', '.7z', '.rar', '.tar.gz')):
        print(f"Skipping zipped folder: {folder}")
        continue
    subj_id = folder.split('_')[0]  # Extract subject ID from folder name
    rest_folder = os.path.join(input_dir, folder, 'MNINonLinear', 'Results', 'rfMRI_REST1_LR')
    target_file = os.path.join(rest_folder, 'rfMRI_REST1_LR_hp2000_clean_rclean_tclean.nii.gz')
    
    if os.path.exists(target_file):
        new_filename = f"{subj_id}_Rest3T.nii.gz"
        os.makedirs(output_dir, exist_ok=True)  # Create target directory if it doesn't exist
        shutil.copy2(target_file, os.path.join(output_dir, new_filename))
        print(f"Copied {target_file} to {os.path.join(output_dir, new_filename)}")
    else:
        print(f"File not found for subject {subj_id}: {target_file}")