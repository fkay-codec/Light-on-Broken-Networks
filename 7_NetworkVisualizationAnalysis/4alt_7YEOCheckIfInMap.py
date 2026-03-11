"""
in this script we import an mni template 

use some RSFC maps

and check if the coordinates from our array are in the RSFC maps
"""
import os
# os.environ["SCIPY_ARRAY_API"] = "1"
from nilearn.image import load_img
import numpy as np
import matplotlib.pyplot as plt
from nilearn.image import coord_transform
from nilearn.plotting import plot_glass_brain, plot_stat_map
from nilearn.datasets import load_mni152_template
from nilearn.image import resample_to_img
import pandas as pd


def extract_dataframes(excel_path, stab_com):
    """
    Input: the excel paths with the module dataset

    Processing: 
        for each gamma create a separate dataframe and check how many points are in the map for each module

    Output: save a text with the results
    """
    df = pd.read_excel(excel_path)
    label = None
    for col in df.columns:
        if f"UC={stab_com}" in col:
            label = col 
            # check how many unique modules there are
            unique_mods=df[col].nunique()
            print(f"\t{col} number of modules: {unique_mods}\n")
    print(label)
    df_new = df[['X', 'Y', 'Z', label, 'ROI']]
    print(df_new.head())

    return df_new

def remap_yeo_labels(network_id):
    """
    remap the yeo numbering to network labeling
    ID	Network name
    1	Visual
    2	Somatomotor
    3	Dorsal Attention
    4	Ventral Attention
    5	Limbic
    6	Frontoparietal
    7	Default

    """
    yeo_network_names = {
        1: "Visual",
        2: "Somatosensory",
        3: "Dorsal Attention",
        4: "Ventral Attention",
        5: "Limbic",
        6: "Frontoparietal",
        7: "Default",
    }
    
    return yeo_network_names.get(network_id, f"Unknown ({network_id})")

def check_how_many_points(dataframes, rsfc_map_path):
    """
    Input: dataframes, SINGLE rsfc_map_path

    Output: txt with the results
    """


    data_sphere=[]
    data_point=[]
    # take the name of the dataframe and store it in new_txt
    label_name = dataframes.columns[3]

    # for each module extract the coordinates
    unique_modules = dataframes[label_name].unique()

    for module in unique_modules:
        df_module = dataframes[dataframes[label_name] == module]

        coordinates = list(zip(df_module['X'], df_module['Y'], df_module['Z']))
        # now we want to check for that module, and these coordinates, how many points are in the rsfc maps

        #! new implementation for single rsfc map
        network_labels_spherical, voxels_in_module = check_sphere_points_multi_network(rsfc_map_path, coordinates=coordinates)
        # network_labels_points
        # Get unique networks and their counts

        # Remap network labels to main network names
        remapped_labels = [remap_yeo_labels(label) for label in network_labels_spherical]

        from collections import Counter
        network_counts = Counter(remapped_labels)

        sum_counts = sum(network_counts.values())

        for network_id, count in sorted(network_counts.items(), key=lambda x: x[1], reverse=True):

            data_sphere.append({
                'Module': module,
                'Network_Name': network_id,
                'Count': count,
                'Prec. of hits': f"{round(100 * (count / sum_counts),2)}%",
                # 'hits/total voxels in module': f"{round(100*(count/voxels_in_module),2)}%"
            })

        network_labels_points = check_points_multi_network(rsfc_map_path, coordinates=coordinates)
        remapped_labels = [remap_yeo_labels(label) for label in network_labels_points]
        network_counts = Counter(remapped_labels)
        sum_counts = sum(network_counts.values())

        for network_id, count in sorted(network_counts.items(), key=lambda x: x[1], reverse=True):
            data_point.append({
                'Module': module,
                'Network_Name': network_id,
                'Count': count,
                'Prec. of hits': f"{round(100 * (count / sum_counts),2)}%",
            })
    out_sphere_df = pd.DataFrame(data_sphere)
    out_point_df = pd.DataFrame(data_point)

    return out_sphere_df, out_point_df
    
def generate_sphere_offsets(radius_voxel):
    """
    Generate a list of voxel offsets that form a sphere of a given radius in voxel space.
    """
    offsets = []
    for x in range(-radius_voxel[0], radius_voxel[0] + 1):
        for y in range(-radius_voxel[1], radius_voxel[1] + 1):
            for z in range(-radius_voxel[2], radius_voxel[2] + 1):
                # Check if the offset is within the sphere radius
                if np.sqrt((x / radius_voxel[0])**2 + (y / radius_voxel[1])**2 + (z / radius_voxel[2])**2) <= 1:
                    offsets.append((x, y, z))
    return offsets
    
def check_sphere_points_multi_network(rsfc_map_path, coordinates=None):
    """
    Input: rsfc_map_path (single map with multiple network labels), coordinates
    Output: list of ALL network labels found across all spheres
    
    Checks ALL network labels within a 5mm sphere around each coordinate.
    Returns ALL networks found (not just most common per sphere).
    """
    if coordinates is None or len(coordinates) == 0:
        raise ValueError("You must provide a list of MNI coordinates.")
    
    # Load the RSFC map
    rsfc_map = load_img(rsfc_map_path)
    glass_brain_template = load_mni152_template()
    rsfc_map_resampled = resample_to_img(rsfc_map, glass_brain_template, interpolation='nearest', copy_header=True, force_resample=True)
    rsfc_data = rsfc_map_resampled.get_fdata()
    
    # Calculate voxel size and sphere offsets
    voxel_size = np.abs(np.diag(rsfc_map_resampled.affine)[:3])
    radius_voxel = (5 / voxel_size).astype(int)
    sphere_offsets = generate_sphere_offsets(radius_voxel)
    
    # Transform coordinates to voxel space
    voxel_coords = []
    for coord in coordinates:
        x_voxel, y_voxel, z_voxel = coord_transform(
            coord[0], coord[1], coord[2], 
            np.linalg.inv(rsfc_map_resampled.affine)
        )
        voxel_coords.append((round(x_voxel), round(y_voxel), round(z_voxel)))

    # Collect ALL network labels from all spheres
    all_network_labels = []
    total_voxels_in_module = 0
    for x_voxel, y_voxel, z_voxel in voxel_coords:
        

        for offset in sphere_offsets:
            x_check = int(x_voxel + offset[0])
            y_check = int(y_voxel + offset[1])
            z_check = int(z_voxel + offset[2])
            
            #Currently we are here: we have voxels both in bounds and out-of-bounds voxels that don't exist in the image

            if (0 <= x_check < rsfc_data.shape[0] and
                0 <= y_check < rsfc_data.shape[1] and
                0 <= z_check < rsfc_data.shape[2]):
                
                total_voxels_in_module +=1 # Counts only valid voxels within image bounds


                label = rsfc_data[x_check, y_check, z_check]
                if label != 0:
                    all_network_labels.append(int(label.item()))
    
    return all_network_labels, total_voxels_in_module

def check_points_multi_network(rsfc_map_path, coordinates=None):
    """
    Input: rsfc_map_path (single map with multiple network labels), coordinates
    Output: list of network labels for each coordinate
    
    Checks which network label each coordinate falls into at the exact MNI point (no sphere).
    """
    if coordinates is None or len(coordinates) == 0:
        raise ValueError("You must provide a list of MNI coordinates.")
    
    # Load the RSFC map
    rsfc_map = load_img(rsfc_map_path)
    glass_brain_template = load_mni152_template()
    rsfc_map_resampled = resample_to_img(rsfc_map, glass_brain_template, interpolation='nearest', copy_header=True, force_resample=True)
    rsfc_data = rsfc_map_resampled.get_fdata()
    
    # Transform coordinates to voxel space
    voxel_coords = []
    for coord in coordinates:
        x_voxel, y_voxel, z_voxel = coord_transform(
            coord[0], coord[1], coord[2], 
            np.linalg.inv(rsfc_map_resampled.affine)
        )
        voxel_coords.append((round(x_voxel), round(y_voxel), round(z_voxel)))
    
    # Check which network each coordinate belongs to (exact point only)
    network_labels = []
    for x_voxel, y_voxel, z_voxel in voxel_coords:
        found_network = 0
        
        # Check only the exact voxel coordinate
        if (0 <= x_voxel < rsfc_data.shape[0] and
            0 <= y_voxel < rsfc_data.shape[1] and
            0 <= z_voxel < rsfc_data.shape[2]):
            
            label = rsfc_data[int(x_voxel), int(y_voxel), int(z_voxel)]
            if label != 0:
                found_network = int(label.item())
        
        network_labels.append(found_network)
    
    return network_labels


stab_coms = 5,6,7

#! implementation using yeo single rsfc map
rsfc_map_yeo_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\yeo\Yeo_JNeurophysiol11_MNI152\Yeo_JNeurophysiol11_MNI152\Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz"

# Get all Excel file paths in the folder
excel_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\7_GraphMetricsGroupLevel_Updated\BrainNetFiles\1_RawDataFrames"
excel_paths = [os.path.join(excel_folder, f) for f in os.listdir(excel_folder) if f.endswith("DF.xlsx")]

for stab_com in stab_coms:
    # for each excel file we get dataframes: 1 dataframe for each gamma that contains: X, Y, Z, Gamma, ROI
    for excel_path in excel_paths:
        print(f"Processing Excel file: {os.path.basename(excel_path)}")
        #! new implementation: extract dataframe with 7 unique communities
        dataframe= extract_dataframes(excel_path, stab_com)
        sphere_df, point_df = check_how_many_points(dataframes=dataframe, rsfc_map_path=rsfc_map_yeo_path)
        excel_name = os.path.basename(excel_path).replace("_Corr_Master_DF.xlsx", "")
        print(f"Saving results for: {excel_name}")
        sphere_out = excel_name + f"_7YEO_Sphere_Check_{stab_com}Com.xlsx"
        point_out = excel_name + f"_7YEO_Point_Check_{stab_com}Com.xlsx"
        # save the dataframes
        directory = os.path.dirname(excel_path)
        sphere_out_path =  os.path.join(directory, sphere_out)
        point_out_path =  os.path.join(directory, point_out)
        sphere_df.to_excel(sphere_out_path, index=False)
        point_df.to_excel(point_out_path, index=False)
        print(f"Results saved to: {sphere_out_path} and {point_out_path}\n")