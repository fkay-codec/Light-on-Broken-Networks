"""
General Function: Plot Heatmap of Correlation Matrix
"""

import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
import numpy as np
sys.path.append(r"C:\Users\foivo\Documents\Coding\Scripts\fNIRS_fMRI_Project")
from utilities import select_folder_path, inverse_z_transform, fisher_z_transform

# this is the template for plotting heatmaps:
def plot_modularity_heatmap(df, title=str, save_path=None, community_labels=None, max=0.8):
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
    # print(df.max().max())
    # quit()
 
    min = -max

    df = df.copy()
    df.columns = [f.replace('_','') for f in df.columns]
    df.index = [f.replace('_','') for f in df.index]
    # print(f"mean value of df: {df.values.mean()}")
    # quit()
    # Spectral_r I dont like it much
    # RdBu_r is better
    # bwr
    # seismic
    # cmap = apply_cmap_alpha('RdBu_r', alpha=0.95)
    # cmap = darken_cmap('bwr', factor=0.8)


    # cmap2 = adjust_gamma_cmap(cmap, gamma=0.9)
    #* this i like: a lot:
    custom = custom_cmap('coolwarm', dark_factor=1.15, alpha=1, gamma=5)

    # custom = custom_cmap('RdBu_r', dark_factor=1, alpha=1, gamma=1)


    plt.figure(figsize=(11, 12))
    heatmap = sns.heatmap(df, annot=False, cmap=custom, center=0, 
                vmin=min, vmax=max, square=True, 
                linewidths=0.05, linecolor='#E6E6FA',

                cbar_kws={'label': "Correlation Coef.", 'shrink': 0.7})
    
    # Move colorbar numbers to the left
    cbar = heatmap.collections[0].colorbar
    cbar.ax.yaxis.tick_left()
    cbar.ax.yaxis.set_label_position('right')
    # cbar.ax.tick_params(length=0)  # Remove the small tick lines

    # Make colorbar tick labels bigger
    cbar.ax.tick_params(labelsize=16)  # Increase font size for tick labels

    cbar.ax.set_ylabel("Correlation Coef.", fontsize=18)  # labelpad controls distance to cbar


    # Move the colorbar position
    pos = cbar.ax.get_position()
    cbar.ax.set_position([pos.x0 + 0.05, pos.y0, pos.width, pos.height])  # Move right by 0.05
    # set the ticks manually
    cbar.set_ticks([min, 0 ,max])# 0.5, 1])  # Set ticks at min, mid, max
    cbar.set_ticklabels([f"{min}<", "0", f"{max}>"])
    
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
    plt.xlabel("Channels", fontsize=12)
    plt.ylabel("Channels", fontsize=12)

    # Set all ticks to show every channel, positioned at center of squares
    plt.xticks([i + 0.5 for i in range(len(df.columns))], df.columns, fontsize=6, rotation=90)
    plt.yticks([i + 0.5 for i in range(len(df.index))], df.index, fontsize=6)
    if not save_path:
        plt.show()

    if save_path:    
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")

max = 0.9
max_name = str(max).replace('.', 'p')

select_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\2_AverageCorM"
out_folder = os.path.join(select_folder, f"Average_Heatmap_{max_name}")
if not os.path.exists(out_folder):
    os.makedirs(out_folder) 
files = [f for f in os.listdir(select_folder) if f.endswith(".xlsx") and "Average" in f]

for file in files:
    input_file = os.path.join(select_folder, file)
    print(f"Processing file: {input_file}")
    df = pd.read_excel(input_file, index_col=0, sheet_name="HbO_Average")
    print(df.head())
    # inverse z-transform for visualization
    df = inverse_z_transform(df)
    print(df.head())



    if "Partial" in file:
        title = "Group Average fNIRS (HbO) FC Matrix \n(Partial Correlations)"
    if "Standard" in file:
        title = "Group Average fNIRS (HbO) FC Matrix \n(Bivariate Correlations)"
    out = out_folder + f"/{file.replace('.xlsx','')}_HbO_Heatmap.png"
    # print(f"Saving heatmap to: {out}")
    # quit()
    plot_modularity_heatmap(df, title=title, save_path=out, max=max)
for file in files:
    input_file = os.path.join(select_folder, file)
    print(f"Processing file: {input_file}")
    df = pd.read_excel(input_file, index_col=0, sheet_name="HbR_Average")
    print(df.head())
    # inverse z-transform for visualization
    df = inverse_z_transform(df)
    print(df.head())



    if "Partial" in file:
        title = "Group Average fNIRS (HbR) FC Matrix \n(Partial Correlations)"
    if "Standard" in file:
        title = "Group Average fNIRS (HbR) FC Matrix \n(Bivariate Correlations)"
    out = out_folder + f"/{file.replace('.xlsx','')}_HbR_Heatmap.png"
    # print(f"Saving heatmap to: {out}")
    # quit()
    plot_modularity_heatmap(df, title=title, save_path=out, max=max)
    