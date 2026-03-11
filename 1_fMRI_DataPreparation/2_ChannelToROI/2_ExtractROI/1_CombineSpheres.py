from nilearn.image import load_img
import nibabel as nib
import numpy as np
from nilearn.masking import apply_mask
from nilearn import plotting, image
import os
import matplotlib.pyplot as plt

# in this script we will iterate through the folder with the binarized 5mm spheres and we want to label them based on their name. ROI1 = label 1, etc.
# then we will combine all the spheres into one nifti file, that containes a different label for each sphere/roi



# Set the directory path
sphere_dir = r"C:\Users\foivo\Desktop\Ubuntu Sharing\Results\1_ROI_Sphere_Creation"  # update this path
sphere_files = [f for f in os.listdir(sphere_dir) if "5mmSphere" in f]
# print(sphere_files)
# quit()
# Initialize combined data array
first_img = load_img(os.path.join(sphere_dir, sphere_files[0]))
combined_data = np.zeros(first_img.shape)


for file in sphere_files:
    # load the sphere
    sphere_img = load_img(os.path.join(sphere_dir, file))

    # extract the ROI number from the filename
    file_name = file.split('_')[0]  # the filename is like: ROI1_S1-D1_5mmSphere.nii.gz || converted to ROI1
    roi_number = int(file_name.replace('ROI', ''))  # extract the number from the filename

    # now label the sphere based on its name
    sphere_data = sphere_img.get_fdata()


    # Check if values are binary if not binarize the data

    # Binarize the sphere data: any value > 0 becomes 1
    sphere_data_binary = (sphere_data > 0).astype(int)
   
    unique_values = np.unique(sphere_data_binary)
    if not (np.array_equal(unique_values, [0., 1.]) or np.array_equal(unique_values, [0]) or np.array_equal(unique_values, [1])):
        print(f"Non-binary values found in {file}: {unique_values}")
        print("Terminating script...")
        quit()

    mask = sphere_data > 0  # create a binary mask
    # Assign the ROI number to the combined data where the mask is True 
    combined_data[mask] = roi_number


# Save the combined labeled image
combined_img = nib.Nifti1Image(combined_data, first_img.affine, first_img.header)
nib.save(combined_img, os.path.join(sphere_dir, "Combined_ROIs.nii.gz"))
print("Combined ROI mask saved with labels.")