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


    plt.figure(figsize=(11, 12))
    heatmap = sns.heatmap(df, annot=False, cmap='coolwarm', center=0, #df.mean().mean(), #df.mean().mean(),
                vmin=vmin, vmax=vmax, square=True, 
                cbar_kws={'label': label, 'shrink': 0.5})
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
                               linewidth=0.4, 
                               alpha=0.5)
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
