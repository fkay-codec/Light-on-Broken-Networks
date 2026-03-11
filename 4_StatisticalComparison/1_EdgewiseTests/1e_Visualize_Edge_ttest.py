import os
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm

# this is the template for plotting heatmaps:
# this is the template for plotting heatmaps:
def plot_ttest_map_categorical(df, title=str, save_path=None, community_labels=None):
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
    # Define p-value boundaries and corresponding colors
    boundaries = [0, 0.001, 0.01, 0.05, 1.0]
    colors = [ "#B2182B", "#EF8A62", "#FDDBC7", "#DADDE8"]  # very red, more red, red, white
    
    # Create discrete colormap
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(boundaries, cmap.N, clip=True)
    
    # Create a mask for the lower triangle (including diagonal)
    mask = np.tril(np.ones_like(df, dtype=bool))

    # plt.figure(figsize=(11, 12))
    # Create the heatmap
    plt.figure(figsize=(11, 12))
    heatmap = sns.heatmap(df, annot=False, cmap=cmap, norm=norm,
                        mask=mask,
                        square=True, 
                        linewidths=0.05, linecolor='#E6E6FA',
                        cbar_kws={'shrink': 0.3})

    # heatmap = sns.heatmap(df, annot=False, cmap=cmap, center=0.05, 
    #             vmin=min, vmax=max, square=True, 
    #             linewidths=0.05, linecolor='#E6E6FA',

    #             cbar_kws={'label': "Correlation Coef.", 'shrink': 0.5})
    
    # Move colorbar numbers to the left
    cbar = heatmap.collections[0].colorbar
    cbar.ax.yaxis.tick_right()
    cbar.ax.yaxis.set_label_position('left')
    # cbar.ax.tick_params(length=0)  # Remove the small tick lines

    # # cbar.ax.tick_params(labelsize=8, length=5, width=2, colors='black')
    
    # cbar.ax.tick_params(length=0, labelsize=8)  # length=0 removes the dashes
    # Remove ALL tick marks
    cbar.ax.tick_params(which='both', length=0, labelsize=12)  # 'which=both' removes major and minor ticks

    # Manually set tick positions at the CENTER of each color band
    cbar.set_ticks([0.0005, 0.0055, 0.03, 0.525])

    # Move the colorbar position
    pos = cbar.ax.get_position()
    cbar.ax.set_position([pos.x0 + 0.025, pos.y0, pos.width, pos.height])  # Move right by 0.05
    # # set the ticks manually
    # cbar.set_ticks([min, 0.05 ,max])# 0.5, 1])  # Set ticks at min, mid, max
    # cbar.set_ticklabels([f"{min}<", "0.05", f"{max}>"])#f"{0.5}", f"{1}"])  # Format tick labels
    # Set custom tick labels
    cbar.set_ticklabels(['p < 0.001', '0.001 ≤ p < 0.01', '0.01 ≤ p < 0.05', 'p ≥ 0.05'], fontsize=12)
    
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
    plt.xticks([i + 0.5 for i in range(len(df.columns))], df.columns, fontsize=4, rotation=90)
    plt.yticks([i + 0.5 for i in range(len(df.index))], df.index, fontsize=4)
 
    if not save_path:
        plt.show()

    if save_path:    
        plt.savefig(save_path, dpi=450, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")
        # quit()


input_fold = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\1_CrossModality_EdgeComparison"
files = [f for f in os.listdir(input_fold) if f.endswith(".xlsx")]

out_folder_name = "Edgewise_ttest_figures_CategCmap"
out_folder = os.path.join(input_fold, out_folder_name)
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

for file in files:
    file_path = os.path.join(input_fold, file)
    clean_name = file.replace(".xlsx", "")

    heme = clean_name.split("_")[-1]
    print(heme)
    if "StandCor" in file:
        title = f"Edgewise Comparison of fMRI vs fNIRS ({heme}):\n(Bivariate Correlations)"
    if "ParCor" in file:
        title = f"Edgewise Comparison of fMRI vs fNIRS ({heme}):\n(Partial Correlations)"
    print(file_path)
    sheet_names = pd.ExcelFile(file_path).sheet_names
    print(sheet_names)
    df = pd.read_excel(file_path, index_col=0, sheet_name=sheet_names[1])
    plot_ttest_map_categorical(df, title, save_path=os.path.join(out_folder, f"{clean_name}_ttest_Cat_map.png"))
    # calculate the % of significant edges
    #extract the upper triangle without the diagonal
    triu_indices = np.triu_indices_from(df, k=1)
    upper_tri_values = df.values[triu_indices]
    sig_edges = np.sum(upper_tri_values < 0.05)
    total_edges = len(upper_tri_values)
    perc_sig = (sig_edges / total_edges) * 100
    print(f"File: {file} - Significant edges: {sig_edges}/{total_edges} ({perc_sig:.2f}%)")
    txt_out_path = os.path.join(out_folder, f"{clean_name}_ttest_significant_edges.txt")
    with open(txt_out_path, "w") as f:
        f.write(f"File: {file}\n")
        f.write(f"Significant edges (p < 0.05): {sig_edges}/{total_edges} ({perc_sig:.2f}%)\n")


