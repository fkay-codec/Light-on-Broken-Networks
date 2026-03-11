import pandas as pd
import numpy as np

def inverse_z_transform(inverse_z_df):
    """Inverse Fisher Z-transformation from the dataframe"""
    # make the diagonal 1
    inversed_z = np.tanh(inverse_z_df)
    np.fill_diagonal(inversed_z.values, 1)
    return inversed_z

def fisher_z_transform(corr_matrix):
    """Apply Fisher z-transformation to a correlation matrix"""
    clip = 0.9999999999  # To avoid infinities
    corr_matrix = corr_matrix.clip(-clip, clip)  # Avoid infinities
    z_matrix = np.arctanh(corr_matrix)
    return z_matrix

# Copy/pasted from Michael Luhrs
from PySide6.QtWidgets import QApplication, QFileDialog
import sys

def select_folder_path(prompt):
    # Open dialog to select the folder of interest and return its path
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)

    # Open a dialog to select a folder 
    folder_path = QFileDialog.getExistingDirectory(None, prompt, r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study")
    return folder_path

def select_file_path(prompt):
    """Open dialog to select the file of interest and return its path."""
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)

    # Open a dialog to select a file
    file_path, _ = QFileDialog.getOpenFileName(None, prompt)
    return file_path

import seaborn as sns
import matplotlib.pyplot as plt
import os
# this is the template for plotting heatmaps:
# this is the template for plotting heatmaps:
def plot_modularity_heatmap(df, title=str, save_path=None, community_labels=None):
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
    max = 0.2
    min = -0.2
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

                cbar_kws={'label': "Correlation Coef.", 'shrink': 0.5})
    
    # Move colorbar numbers to the left
    cbar = heatmap.collections[0].colorbar
    cbar.ax.yaxis.tick_left()
    cbar.ax.yaxis.set_label_position('left')
    # cbar.ax.tick_params(length=0)  # Remove the small tick lines

    # Move the colorbar position
    pos = cbar.ax.get_position()
    cbar.ax.set_position([pos.x0 + 0.05, pos.y0, pos.width, pos.height])  # Move right by 0.05
    # set the ticks manually
    cbar.set_ticks([min, 0 ,max])# 0.5, 1])  # Set ticks at min, mid, max
    cbar.set_ticklabels([f"{min}<", "0", f"{max}>"])#f"{0.5}", f"{1}"])  # Format tick labels
    
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
    plt.xlabel("Channels", fontsize=8)
    plt.ylabel("Channels", fontsize=8)

    # Set all ticks to show every channel, positioned at center of squares
    plt.xticks([i + 0.5 for i in range(len(df.columns))], df.columns, fontsize=4, rotation=90)
    plt.yticks([i + 0.5 for i in range(len(df.index))], df.index, fontsize=4)
    if not save_path:
        plt.show()

    if save_path:    
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")
     
