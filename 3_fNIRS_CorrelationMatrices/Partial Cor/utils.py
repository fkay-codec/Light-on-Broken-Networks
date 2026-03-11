import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_heatmap(df, title=str, save_path=None, community_labels=None):
    """
    Plot a heatmap for a 2D DataFrame (typically correlation matrices).
    
    Parameters:
    -----------
    df : pandas.DataFrame
        2D DataFrame to plot as heatmap
    title : str
        Title for the plot
    """
    if "Correlation" in title:
        actual_min = df.min().min()
        actual_max = df.max().max()
        vmin, vmax = actual_min, actual_max
        label = 'Correlation Coefficient'
    elif "T-test" in title or "t-statistic" in title.lower():
        # For t-statistics, use data range but center at 0
        data_max = max(abs(df.min().min()), abs(df.max().max()))
        vmin, vmax = -data_max, data_max
        label = 'T-statistic'
    else:
        vmin, vmax = None, None
        label = 'Value'

    # Create mask for diagonal to make it black
    mask = np.eye(df.shape[0], dtype=bool)

    plt.figure(figsize=(11, 12))
    heatmap = sns.heatmap(df, annot=False, cmap='RdBu_r', center=None,#df.mean().mean(), #df.mean().mean(), #df.mean().mean(),
                vmin=vmin, vmax=vmax, square=True, 
                cbar_kws={'label': label, 'shrink': 0.5},
                mask=mask,  # Apply the mask to hide diagonal
    )
    # Manually draw black squares on the diagonal
    for i in range(df.shape[0]):
        heatmap.add_patch(plt.Rectangle((i, i), 1, 1, fill=True, color='black'))

    # Move colorbar numbers to the left
    cbar = heatmap.collections[0].colorbar
    cbar.ax.yaxis.tick_left()
    cbar.ax.yaxis.set_label_position('left')
    # cbar.ax.tick_params(length=0)  # Remove the small tick lines

    # Move the colorbar position
    pos = cbar.ax.get_position()
    cbar.ax.set_position([pos.x0 + 0.05, pos.y0, pos.width, pos.height])  # Move right by 0.05
    cbar.set_ticks([vmin, (vmin + vmax) / 2, vmax])  # Set ticks at min, mid, max
    cbar.set_ticklabels([f"{vmin:.2f}", f"{(vmin + vmax) / 2:.2f}", f"{vmax:.2f}"])  # Format tick labels
    # # Reduce number of colorbar ticks
    # if "Correlation" in title:
    #     cbar.set_ticks([-1, -0.5, 0, 0.5, 1])  # 5 ticks for correlation
    # elif "T-test" in title or "t-statistic" in title.lower():
    #     data_max = max(abs(df.min().min()), abs(df.max().max()))
    #     cbar.set_ticks([-data_max, -data_max/2, 0, data_max/2, data_max])  # 5 ticks for t-stats
    # else:
    #     cbar.locator = plt.MaxNLocator(5)  # Default to 5 ticks for other cases
  
    # Add community boundary lines if community labels are provided

    ### This creates a chessboard pattern which i dont like
    # if community_labels is not None:
    #     # Find boundaries between communities
    #     boundaries = []
    #     current_community = community_labels[0]
        
    #     for i, community in enumerate(community_labels[1:], 1):
    #         if community != current_community:
    #             boundaries.append(i)
    #             current_community = community
        
    #     # Draw boundary lines ONLY at community transitions
    #     for boundary in boundaries:
    #         # Vertical lines
    #         plt.axvline(x=boundary, color='grey', linewidth=0.4, alpha=0.5)
    #         # Horizontal lines
    #         plt.axhline(y=boundary, color='grey', linewidth=0.4, alpha=0.5)
        
    #     print(f"Added {len(boundaries)} community boundary lines at positions: {boundaries}")
    
    ### this creates rectangles around each community
    if community_labels is not None:
        # Find boundaries between communities
        boundaries = [0]  # Start at 0
        current_community = community_labels[0]
        
        for i, community in enumerate(community_labels[1:], 1):
            if community != current_community:
                boundaries.append(i)
                current_community = community
        
        boundaries.append(len(community_labels))  # End boundary
        
        print(f"Community boundaries at positions: {boundaries}")
        
        # Draw rectangles around each community
        for i in range(len(boundaries) - 1):
            start = boundaries[i]
            end = boundaries[i + 1]
            
            # Draw rectangle around this community
            # Rectangle coordinates: (x_start, y_start, width, height)
            rect = plt.Rectangle((start, start), 
                               end - start,  # width
                               end - start,  # height
                               fill=False,  # No fill, just outline
                               edgecolor='grey', 
                               linewidth=1.5, 
                               alpha=0.8)
            plt.gca().add_patch(rect)
            
            print(f"Added rectangle for community from {start} to {end}")


    plt.title(title)
    plt.xlabel("Channels", fontsize=8)
    plt.ylabel("Channels", fontsize=8)

    # Set all ticks to show every channel, positioned at center of squares
    plt.xticks([i + 0.5 for i in range(len(df.columns))], df.columns, fontsize=4, rotation=90)
    plt.yticks([i + 0.5 for i in range(len(df.index))], df.index, fontsize=4)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")
     
    plt.show()

def plot_correlation_2d_heatmap(df, title=str, save_path=None):
    """Convenience function specifically for correlation matrices."""
    plot_heatmap(df, title=title, save_path=save_path)

def plot_2d_heatmap(df, title=str, save_path=None):
    """Convenience function for any 2D DataFrame."""
    plot_heatmap(df, title=title, save_path=save_path)

### Folder/File Selecting scripts
# Copy/pasted from Michael Luhrs

from PySide6.QtWidgets import QApplication, QFileDialog
import sys

def select_folder_path(prompt):
    # Open dialog to select the folder of interest and return its path
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)

    # Open a dialog to select a folder 
    folder_path = QFileDialog.getExistingDirectory(None, prompt, r"C:\Users\foivo\Documents\Satori\SampleData\Tilburg_Paper")
    return folder_path

def select_file_path(prompt):
    """Open dialog to select the file of interest and return its path."""
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)

    # Open a dialog to select a file
    file_path, _ = QFileDialog.getOpenFileName(None, prompt, r"C:\Users\foivo\Documents\Satori\SampleData\Tilburg_Paper")
                                               
    return file_path


if __name__ == "__main__":
# simple test of the heatmap function with small and large matrices

    # Create a test 2D correlation matrix
    columns = ['a', 'b', 'c', 'd', 'e']
    data = [[1.0, 0.8, -0.2, 0.1, 0.4],
            [0, 1.0, -0.1, 0.3, 0.2],
            [0, -0.1, 1.0, 0.6, -0.5],
            [0, 0.3, 0.6, 1.0, 0.7],
            [0, 0.2, -0.5, 0.7, 1.0]]

    test_df = pd.DataFrame(data, columns=columns, index=columns)
    print("Test correlation matrix:")
    print(test_df)
    # quit()


    # Test the heatmap function
    plot_heatmap(test_df, "Test Correlation Matrix")

    # create a test 2D cor matrix with 100 channels
    columns_100 = [f'ch_{i+1}' for i in range(100)]
    data_100 = np.random.rand(100, 100)*2 - 1  # random values between -1 and 1
    data_100 = (data_100 + data_100.T) / 2  # make it symmetric
    np.fill_diagonal(data_100, 1)  # set diagonal to 1
    test_df_100 = pd.DataFrame(data_100, columns=columns_100, index=columns_100)
    print("Test correlation matrix with 100 channels:")
    print(test_df_100)
    # Test the heatmap function with 100 channels
    plot_heatmap(test_df_100, "Test Correlation Matrix with 100 Channels")