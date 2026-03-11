
"""
NOTE: Check with util plot_heat_function to see if there is something u need to incorporate
In this script we are going to sort the matrices based on their modularity.
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append(r"C:\Users\foivo\Documents\Coding\Scripts\fNIRS_fMRI_Project")
import utilities as util

# def mytitle(title):
#     """
#     Change the title to your liking
#     """
#     if "partial" in title.lower():
#         title = ""
#     elif "standard" in title.lower():
#         title = ""
#     return title
def inverse_z_transform(inverse_z_df):
    """Inverse Fisher Z-transformation from the dataframe"""
    # make the diagonal 1
    inversed_z = np.tanh(inverse_z_df)
    np.fill_diagonal(inversed_z.values, 1)
    return inversed_z

def sort_matrix_by_modularity(matrix, labels):
    """
    Input:
    - matrix: 2D Dataframe with correlation matrix 
    - labels: 1D pandas Series with the labels of each node
    Output:
    - sorted_matrix: 2D Dataframe with the sorted correlation matrix
    """
    # find the number of unique labels
    unique_labels, unique_counts = np.unique(labels, return_counts=True) 

    # sort them by their size
    # get the indices that would sort them by size [largest first]
    idx = np.argsort(unique_counts)[::-1]

    # now reorder the Pandas Series based on the sorted indices
    sorted_labels = unique_labels[idx]

    # Create a mapping of the original labels to their new sorted order
    label_order_mapping = {label: i for i, label in enumerate(sorted_labels)}

    # Map the original labels to their new order
    labels_sorted = labels.map(label_order_mapping)
    labels_sorted = labels_sorted.sort_values()
    # now that we have the labels sorted by size and community, we want to sort each community by their name so they have a specific order every time

    # create a dataframe with the sorted labels, and their names on a new column
    df = pd.DataFrame({
        labels_sorted.name: labels_sorted,
        'ROI': labels_sorted.index
    })
    # manipulate the channel column to remove the stuff after the underscore
    df['ROI'] = df['ROI'].str.split('_').str[0]
    df['ROI'] = df['ROI'].str.replace('ROI', '').astype(int)
    # now for each community, sort the ROIs by their number
    df= df.sort_values(by=[labels_sorted.name, 'ROI'])
    final_labels_series = df.iloc[:,0]
    # print(final_labels_series)
    # quit()
    # now we get this index as the new sorted labels
    # Sort the rows and columns of the matrix based on the label order
    sorted_indices = df.index

    sorted_matrix = matrix.loc[sorted_indices, sorted_indices]
    return sorted_matrix, final_labels_series

def plot_modularity_heatmap(df, title=str, save_path=None, community_labels=None, max=0.8, uc=None):
    from matplotlib.colors import ListedColormap
    import matplotlib as mpl
    def apply_cmap_alpha(cmap_name='RdBu_r', alpha=0.85):
        base_cmap = mpl.colormaps.get_cmap(cmap_name)
        cmap_colors = base_cmap(np.arange(base_cmap.N))
        cmap_colors[:, -1] = alpha  # modify alpha channel
        return ListedColormap(cmap_colors)  
    def darken_cmap(cmap_name='bwr', factor=0.8):
        base = mpl.colormaps.get_cmap(cmap_name)
        colors = base(np.linspace(0, 1, base.N))
        colors[:, :3] = np.clip(colors[:, :3] * factor, 0, 1)
        return ListedColormap(colors)
    def adjust_gamma_cmap(cmap_name='bwr', gamma=0.8):
        base = mpl.colormaps.get_cmap(cmap_name)
        colors = base(np.linspace(0, 1, base.N))
        colors[:, :3] = np.power(colors[:, :3], gamma)
        return ListedColormap(colors)
    
    def custom_cmap(cmap_name='bwr', dark_factor=0.8, alpha=0.9, gamma=0.9):
        base = mpl.colormaps.get_cmap(cmap_name)
        colors = base(np.linspace(0, 1, base.N))
        colors[:, :3] = np.clip(colors[:, :3] * dark_factor, 0, 1)
        colors[:, :3] = np.power(colors[:, :3], gamma)
        colors[:, -1] = alpha
        return ListedColormap(colors)

        
    """
    Plot a heatmap for a 2D DataFrame (typically correlation matrices).
    
    Parameters:
    -----------
    df : pandas.DataFrame
        2D DataFrame to plot as heatmap
    title : str
        Title for the plot
    """

    
    min = -max

    # cmap2 = adjust_gamma_cmap(cmap, gamma=0.9)
    #* this i like: a lot:
    custom = custom_cmap('coolwarm', dark_factor=1.15, alpha=1, gamma=5)

    # custom = custom_cmap('RdBu_r', dark_factor=1, alpha=1, gamma=1)

    df = df.copy()
    df.columns = [f.replace('_','/') for f in df.columns]
    df.index = [f.replace('_','/') for f in df.index]

    plt.figure(figsize=(11, 12))
    heatmap = sns.heatmap(df, annot=False, cmap=custom, center=0, 
                vmin=min, vmax=max, square=True, 
                linewidths=0.05, linecolor='#E6E6FA',

                cbar_kws={'label': "Normalized Correlation Coef.", 'shrink': 0.6})
    
    # Move colorbar numbers to the left
    cbar = heatmap.collections[0].colorbar
    cbar.ax.yaxis.tick_left()
    cbar.ax.yaxis.set_label_position('right')
    # cbar.ax.tick_params(length=0)  # Remove the small tick lines

    # Make colorbar tick labels bigger
    cbar.ax.tick_params(labelsize=12)  # Increase font size for tick labels

    cbar.ax.set_ylabel("Correlation Coef. (Normalized)", fontsize=18)  # labelpad controls distance to cbar


    # Move the colorbar position
    pos = cbar.ax.get_position()
    cbar.ax.set_position([pos.x0 + 0.1, pos.y0, pos.width, pos.height])  # Move right by 0.1
    # set the ticks manually
    cbar.set_ticks([min, 0 ,max])# 0.5, 1])  # Set ticks at min, mid, max
    cbar.set_ticklabels([f"{min}", "0", f"{max}"])

    
    ### this creates rectangles around each community
    if community_labels is not None:
        # Find boundaries between communities
        boundaries = [0]  # Start at 0
        current_community = community_labels[0]
        print(community_labels)
        print(community_labels.name)
        print(current_community)
        # quit()
        for i, community in enumerate(community_labels):
            if community != current_community:
                boundaries.append(i)
                boundarie_name = community_labels.index[i]
                print(f"Community changed at index {i}, row name {boundarie_name}")
                current_community = community
            print(f"Index: {i}, Community: {community}, Current: {current_community} with row name {community_labels.index[i]}")
        print(f"Final boundaries: {boundaries} with length {len(boundaries)}")

        boundaries.append(len(community_labels))  # End boundary
        
        print(f"Community boundaries at positions: {boundaries}")
        # Draw rectangles around each community
        for i in range(len(boundaries)-1):
            start = boundaries[i]
            end = boundaries[i + 1]

            
            # Draw rectangle around this community
            # Rectangle coordinates: (x_start, y_start, width, height)
            rect = plt.Rectangle((start, start), 
                               end - start,  # width
                               end - start,  # height
                               fill=False,  # No fill, just outline
                               edgecolor='black', 
                               linewidth=2.3, 
                               alpha=0.5)
            plt.gca().add_patch(rect)
            
            print(f"Added rectangle for community from {start} to {end}")
            # quit()

    plt.title(title, fontsize=20, fontweight='bold')
    if "fNIRS" in title:
        m_label = "Channels"
    elif "fMRI" in title:
        m_label = "ROIs"
    plt.xlabel(m_label, fontsize=12)
    plt.ylabel(m_label, fontsize=12)


    # Set all ticks to show every channel, positioned at center of squares
    plt.xticks([i + 0.5 for i in range(len(df.columns))], df.columns, fontsize=6, rotation=90)
    plt.yticks([i + 0.5 for i in range(len(df.index))], df.index, fontsize=6)

    plt.tight_layout()
    if save_path:
        sanitized_title = title.replace(" = ","").replace(" ","_").replace(".","p").replace("\n","").replace("(","").replace(")","").replace("=",'').replace(";","").replace(":","_")

        if uc is not None:
            out_filename =os.path.join(save_path, sanitized_title + f"_UC{uc}.png")
        else:
            out_filename =os.path.join(save_path, sanitized_title + ".png")

        plt.savefig(out_filename, dpi=450, bbox_inches='tight')
        print(f"Figure saved to: {out_filename}")
     


def process_fmri_files(data_path, label_path, name='', out_save = None, max=0.8):
    data = pd.read_excel(data_path, index_col=0)
    #! the input data here where first z-transformed. then normalized across matrices HbO/HbR/fMRI to be in the same linear scale.
    #! therefore there is no need to invert them back to correlation values r because they are not in correlation scale anymore but rather (Normalized) r
    
    # # invert them back to correlation values r
    # data = inverse_z_transform(data)

    labels = pd.read_excel(label_path, index_col=0, sheet_name='Modularity')

    column_names = labels.columns.tolist()
    for col in range(len(column_names)):

        current_column = labels[column_names[col]] # we get the first column as a pandas series

        sorted_df, labels_sorted = sort_matrix_by_modularity(data, current_column)
        col_name = current_column.name

        q, gamma, uc = col_name.split('|')
        gamma = float(gamma.replace('G=', ''))
        q = float(q.replace('Q=', ''))
        uc = int(uc.replace('UC=', ''))
        print(f"Gamma: {gamma}, Q: {q}")
        title = f"fMRI FC Modularity Matrix\n({name})"


        # quit()
        plot_modularity_heatmap(sorted_df,
                                title=title,
                                save_path=out_save,
                                community_labels=labels_sorted,
                                max=max,
                                uc = uc
        )

def rename_to_standard(df_original_list):
    """We want the output correlation matrices to have matching labeling system therefore we are creating a universal function for bothf NIRS and MRI data for this project"""
    reference_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\1_ROI_Sphere_Creation"
    # this path contains the original mapping scheme which exists in the nametags of the files
    file_names_in_reference = [f for f in os.listdir(reference_path) if 'point' in f and f.startswith('ROI')]
    cleaned_file_names = [name.split('_point')[0].replace('-','') for name in file_names_in_reference]
    roi_list = [f.split("_")[0] for f in cleaned_file_names]
    channel_list = [f.split("_")[1] for f in cleaned_file_names]

    # replace underscores
    df_channel_list = df_original_list.copy()
    df_channel_list = [ch.replace('_','') for ch in df_channel_list]

    # check if it is indeed the right mapping
    for roi, ch, original in zip(roi_list, channel_list, cleaned_file_names):
        constructed_name = f"{roi}_{ch}"
        if constructed_name != original:
            print(f"Mismatch found: constructed {constructed_name} vs original {original}") 
            quit()
    # now we want to create this mapping from roi_list = channel_list = cleaned_file_names
    roi_channel_map = dict(zip(roi_list, channel_list))
    for roi, channel in roi_channel_map.items():
        # print(f"Mapping ROI: {roi} to Channel: {channel}")
        # quit()
        for df_ch in df_channel_list:
            # print(f"Comparing {channel} with {df_ch}")
            if channel == df_ch:
                # print(f"found that {roi}, {channel} matches {df_ch}")
                # reconstruct its name now
                new_name = f"{roi}_{channel}"
                # print(f"Renaming {df_ch} to {new_name}")
                df_channel_list[df_channel_list.index(df_ch)] = new_name
                
                # quit()
    # before we return we want to check if it is indeed correct, the ordering 

    verify = df_original_list.copy()
    verify = [ch.replace('_','') for ch in verify]
    verify_with = df_channel_list.copy()
    verify_with = [ch.split('_')[1] for ch in verify_with]
    for v1, v2 in zip(verify, verify_with):
        if v1 != v2:
            print(f"Final verification mismatch: {v1} vs {v2}")
            quit()
    print("Renaming verification successful.")
    # quit()

    return df_channel_list

def process_fnirs_files(data_path, label_path, name='', out_save = None, max=0.8):
    dataframe = pd.ExcelFile(data_path)
    dataframe_sheet_names = dataframe.sheet_names
    print(f"Data sheets found: {dataframe_sheet_names}")
    # quit()
    # we want to do it for both data sheets
    for sheet in dataframe_sheet_names:
        data = pd.read_excel(dataframe, sheet_name=sheet, index_col=0)
        channels_for_data = rename_to_standard(data.index.tolist())

        data.columns = channels_for_data
        data.index = channels_for_data

        print(f"Processing sheet: {sheet}")
        heme = sheet.split('_')[0]
        # invert them back to correlation values r
        # data = inverse_z_transform(data)

        # find the corresponding label sheet
        labels_dataframe = pd.ExcelFile(label_path)
        label_sheet_names = labels_dataframe.sheet_names
        for label_sheet in label_sheet_names:
            if sheet.split('_')[0] in label_sheet and 'Modularity' in label_sheet:
                print(f"Found matching label sheet: {label_sheet} for data sheet: {sheet}")
                correct_label_sheet = label_sheet
                break
        labels = pd.read_excel(label_path, index_col=0, sheet_name=correct_label_sheet)
 
        channels_for_labels = rename_to_standard(labels.index.tolist())

        # extra verification for me only:
        verify = labels.index.tolist()
        verify = [ch.replace('_','') for ch in verify]
        verify_with = channels_for_data.copy()
        verify_with = [ch.split('_')[1] for ch in verify_with]
        for v1, v2 in zip(verify, verify_with):
            if v1 != v2:
                print(f"Final verification mismatch: {v1} vs {v2}")
                quit()
        print("Renaming verification successful.")

        extra_verify = channels_for_data.copy()
        extra_verify_with = channels_for_labels.copy()
        for v1, v2 in zip(extra_verify, extra_verify_with):
            if v1 != v2:
                print(f"Channel mismatch: {v1} vs {v2}")
                quit()
        labels.index = channels_for_labels


        column_names = labels.columns.tolist()
        for col in range(len(column_names)):

            current_column = labels[column_names[col]] # we get the first column as a pandas series

            sorted_df, labels_sorted = sort_matrix_by_modularity(data, current_column)
            col_name = current_column.name

            q, gamma, uc = col_name.split('|')
            uc = int(uc.replace('UC=', ''))
            gamma = float(gamma.replace('G=', ''))
            q = float(q.replace('Q=', ''))
            print(f"Gamma: {gamma}, Q: {q}")
            title = f"fNIRS ({heme}) FC Modularity Matrix\n({name})"
            # quit()
            plot_modularity_heatmap(sorted_df,
                                    title=title,
                                    save_path=out_save,
                                    community_labels=labels_sorted,
                                    max=max,
                                    uc= uc
            )

def process_fnirs_files_with_fmri_label(data_path, label_path, name='', out_save = None, max=0.8):
    dataframe = pd.ExcelFile(data_path)
    dataframe_sheet_names = dataframe.sheet_names
    print(f"Data sheets found: {dataframe_sheet_names}")
    # quit()
    # we want to do it for both data sheets
    for sheet in dataframe_sheet_names:
        data = pd.read_excel(dataframe, sheet_name=sheet, index_col=0)

        channels = data.index.tolist()
        channels_for_data = rename_to_standard(channels)
        data.columns = channels_for_data
        data.index = channels_for_data

        print(f"Processing sheet: {sheet}")
        heme = sheet.split('_')[0]
        # invert them back to correlation values r
        # data = inverse_z_transform(data)
        # find the corresponding label sheet

        labels = pd.read_excel(label_path, index_col=0, sheet_name='Modularity')
        channels_for_labels = labels.index.tolist()
        print(channels_for_data)
        print(channels_for_labels)
        # before the check sort the channles for data and channels for labels
        # HARD CHECK

        if len(channels_for_data) != len(channels_for_labels):
            print(f"Channel length mismatch: {len(channels_for_data)} vs {len(channels_for_labels)}")
            quit()
        for ch1, ch2 in zip(channels_for_data, channels_for_labels):
            if ch1 != ch2:
                print(f"Channel mismatch: {ch1} vs {ch2}")
                quit()

        labels.index = channels_for_labels

        column_names = labels.columns.tolist()
        for col in range(len(column_names)):

            current_column = labels[column_names[col]] # we get the first column as a pandas series

            sorted_df, labels_sorted = sort_matrix_by_modularity(data, current_column)
            col_name = current_column.name

            q, gamma, uc= col_name.split('|')
            uc = int(uc.replace('UC=', ''))
            gamma = float(gamma.replace('G=', ''))
            q = float(q.replace('Q=', ''))
            print(f"Gamma: {gamma}, Q: {q}")
            title = f"Group Average fNIRS ({heme}) FC Matrix Sorted by fMRI Modularity\n({name})"
            # quit()
            plot_modularity_heatmap(sorted_df,
                                    title=title,
                                    save_path=out_save,
                                    community_labels=labels_sorted,
                                    max=max,
                                    uc= uc
            )


# max = 0.9
max =1.0
max_name = str(max).replace('.','p')

out_folder_name = f"Modularity_Figures_{max_name}"
#! use the normalized or the other for the data? If we use the normalized then the scaling will be the same between the two and we wont have to worry about discrepancies in terms of r-value ranges
fmri_data_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\5_AverCorM_NormForGraph" 
fmri_label_folder =r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\7_GraphMetricsGroupLevel_Updated"

#! # alternative folder
# fmri_label_folder=r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\7_GraphMetricsGroupLevelMaxQ"


fmri_data_files = [f for f in os.listdir(fmri_data_folder) if f.endswith('.xlsx') and 'Average' in f]
fmri_label_files = [f for f in os.listdir(fmri_label_folder) if f.endswith('.xlsx') and 'GraphMetrics' in f]

for data, label in zip(fmri_data_files, fmri_label_files):
    if 'StandardCor' in data and 'StandardCor' in label:
        data_path = os.path.join(fmri_data_folder, data)
        label_path = os.path.join(fmri_label_folder, label)
        print(f"Processing data file: {data_path}")
        print(f"With label file: {label_path}")
        out_save = os.path.join(fmri_label_folder, out_folder_name)
        if not os.path.exists(out_save):
            os.makedirs(out_save)
        process_fmri_files(data_path, label_path, name='Bivariate Correlations', out_save = out_save, max=max)
    if 'PartialCor' in data and 'PartialCor' in label:
        data_path = os.path.join(fmri_data_folder, data)
        label_path = os.path.join(fmri_label_folder, label)
        out_save = os.path.join(fmri_label_folder, out_folder_name)
        if not os.path.exists(out_save):
            os.makedirs(out_save)
        print(f"Processing data file: {data_path}")
        print(f"With label file: {label_path}")
        process_fmri_files(data_path, label_path, name='Partial Correlations',out_save = out_save, max=max)


# For fNIRS data we are going to do 2 things:
# 1) process the data and sort them based on fNIRS modularity
# 2) process the data and sort them based on fMRI modularity and save them in a different folder alltogether

fnirs_data_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\3_AverCorM_NormForGraph"
fnirs_label_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\5_GraphMetricsGroupLevel_Updated"

#! # alternative folder
# fnirs_label_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\5_GraphMetricsGroupLevelMaxQ"

fnirs_data_files = [f for f in os.listdir(fnirs_data_folder) if f.endswith('.xlsx') and 'Average' in f]
fnirs_label_files = [f for f in os.listdir(fnirs_label_folder) if f.endswith('.xlsx') and 'GraphMetrics' in f]

out_folder_name_for_fmri_sorted = f"Modularity_Figures_Sorted_BasedOnFMRI"

for data, label in zip(fnirs_data_files, fnirs_label_files):
    if 'StandardCor' in data and 'StandardCor' in label:
        data_path = os.path.join(fnirs_data_folder, data)
        label_path = os.path.join(fnirs_label_folder, label)
        out_save = os.path.join(fnirs_label_folder, out_folder_name)
        if not os.path.exists(out_save):
            os.makedirs(out_save)
        print(f"Processing data file: {data_path}")
        print(f"With label file: {label_path}")
        process_fnirs_files(data_path, label_path, name='Bivariate Correlations', out_save = out_save, max=max)
    if 'PartialCor' in data and 'PartialCor' in label:
        data_path = os.path.join(fnirs_data_folder, data)
        label_path = os.path.join(fnirs_label_folder, label)
        out_save = os.path.join(fnirs_label_folder, out_folder_name)
        if not os.path.exists(out_save):
            os.makedirs(out_save)
        print(f"Processing data file: {data_path}")
        print(f"With label file: {label_path}")
        process_fnirs_files(data_path, label_path, name='Partial Correlations', out_save = out_save, max=max)




# 2) now we are going to sort the fNIRS data based on the fMRI modularity and save them to the below folder
out_folder_name_for_fmri_sorted = f"Modularity_Figures_Sorted_BasedOnFMRI"

for data, label in zip(fnirs_data_files, fmri_label_files):
    if 'StandardCor' in data and 'StandardCor' in label:
        print(f"Processing fNIRS data file: {data} with fMRI label file: {label}")
        data_path = os.path.join(fnirs_data_folder, data)
        label_path = os.path.join(fmri_label_folder, label)
        print(f"with their paths: \nData:{data_path}\nLabel:{label_path}")
        out_save = os.path.join(fnirs_label_folder, out_folder_name_for_fmri_sorted)
        if not os.path.exists(out_save):
            os.makedirs(out_save)
        process_fnirs_files_with_fmri_label(data_path, label_path, name='Bivariate Correlations', out_save = out_save, max=max)
    if "PartialCor" in data and 'PartialCor' in label:
        print(f"Processing fNIRS data file: {data} with fMRI label file: {label}")
        data_path = os.path.join(fnirs_data_folder, data)
        label_path = os.path.join(fmri_label_folder, label)
        print(f"with their paths: \nData:{data_path}\nLabel:{label_path}")
        out_save = os.path.join(fnirs_label_folder, out_folder_name_for_fmri_sorted)
        if not os.path.exists(out_save):
            os.makedirs(out_save)
        process_fnirs_files_with_fmri_label(data_path, label_path, name='Partial Correlations', out_save = out_save, max=max)