"""
0.55 0.0 0.0	% Vivid Red
0.10 0.70 0.20    % Emerald Green
0.20 0.40 0.90    % Royal Blue
0.98 0.85 0.20    % Gold Yellow
0.80 0.20 0.80    % Violet
0.10 0.80 0.80    % Turquoise
1.00 0.50 0.00    % Orange
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerPatch

# Set Times New Roman as default font
plt.rcParams['font.family'] = 'Times New Roman'

class HandlerCircle(HandlerPatch):
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = mpatches.Circle(xy=center, radius=height/2, 
                          facecolor=orig_handle.get_facecolor(),
                          edgecolor=orig_handle.get_edgecolor(),
                          linewidth=orig_handle.get_linewidth())
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]

def create_legend(colors, title, output_path_base, with_border=True):
    """Create a legend with specified colors and save it."""
    fig, ax = plt.subplots(figsize=(12, 14), dpi=300)
    
    # Create legend handles with circles
    legend_handles = [mpatches.Circle((0, 0), radius=1, facecolor=color, edgecolor='black', linewidth=1) 
                      for color, label in colors]
    
    # Create legend with custom labels
    legend_labels = [label for color, label in colors]
    
    # Add legend to the plot
    legend = ax.legend(legend_handles, legend_labels, 
              loc='center', 
              frameon=with_border, 
              fontsize=56,
              title=title,
              title_fontsize=56,
              ncol=1,
              handlelength=2.5,
              handleheight=2.5,
              handler_map={mpatches.Circle: HandlerCircle()},
              prop={'weight': 'bold', 'size': 56})
    
    # Center the legend title and make it bold
    legend.get_title().set_ha('center')
    legend.get_title().set_weight('bold')
    
    # Remove axes
    ax.axis('off')
    
    # Adjust layout
    # plt.tight_layout()
    
    # Save with appropriate suffix
    border_suffix = "_border" if with_border else "_no_border"
    output_path = output_path_base.replace('.png', f'{border_suffix}.png')
    plt.savefig(output_path, dpi=300, bbox_inches=None, pad_inches=0)
    print(f"Figure saved to: {output_path}")
    
    plt.close()


#* for fMRI bivariate legend

fmri_biv_colors = [
    ([0.55, 0.0, 0.0], r'1$^{st}$: DMN'),
    ([0.10, 0.70, 0.20], r'2$^{nd}$: CEN, SN'),
    ([0.20, 0.40, 0.90], r'3$^{rd}$: d-AN'),
    ([0.98, 0.85, 0.20], r'4$^{th}$: d-AN'),
    ([0.80, 0.20, 0.80], r'5$^{th}$: SMN'),
    ([0.10, 0.80, 0.80], r'6$^{th}$: p-VN, c-VN'),
    ([1.00, 0.50, 0.00], r'7$^{th}$: c-VN')
]

fmri_biv_output_base = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\99_figures_Gimp\BrainNet Figz\fMRI_Biv_ColorMap.png"

# Create versions with and without border
# create_legend(fmri_biv_colors, "fMRI Subnetworks", fmri_biv_output_base, with_border=True)
create_legend(fmri_biv_colors, "fMRI Subnetworks", fmri_biv_output_base, with_border=False)


#* for fNIRS HbO/HbR Bivariate

fnirs_hbo_colors = [
    ([0.55, 0.0, 0.0], r'1$^{st}$'), # vivid red
    ([0.10, 0.70, 0.20], r'2$^{nd}$'), # emerald green
    # ([0.20, 0.40, 0.90], r'3$^{rd}$'), # royal blue
    ([0.80, 0.20, 0.80], r'4$^{th}$'), # violet
    ([0.98, 0.85, 0.20], r'5$^{th}$'), # gold yellow
    ([0.20, 0.40, 0.90], r'6$^{th}$'), # royal blue

    # ([0.10, 0.80, 0.80], r'6$^{th}$'), # turquoise
    ([1.00, 0.50, 0.00], r'7$^{th}$') # orange
]

fnirs_hbo_output_base = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\99_figures_Gimp\BrainNet Figz\fNIRS_Biv_HbO_ColorMap.png"

# Create versions with and without border
# create_legend(fnirs_hbo_colors, "Subnetworks", fnirs_hbo_output_base, with_border=True)
create_legend(fnirs_hbo_colors, "Subnetworks", fnirs_hbo_output_base, with_border=False)

#hbr
fnirs_hbr_colors = [
    ([0.55, 0.0, 0.0], r'1$^{st}$'), # vivid red
    ([0.10, 0.70, 0.20], r'2$^{nd}$'), # emerald green
    # ([0.20, 0.40, 0.90], r'3$^{rd}$'), # royal blue
    ([0.80, 0.20, 0.80], r'4$^{th}$'), # violet
    # ([0.98, 0.85, 0.20], r'5$^{th}$'), # gold yellow
    ([0.20, 0.40, 0.90], r'6$^{th}$'), # royal blue

    # ([0.10, 0.80, 0.80], r'6$^{th}$'), # turquoise
    ([1.00, 0.50, 0.00], r'7$^{th}$') # orange
]

fnirs_hbr_output_base = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\99_figures_Gimp\BrainNet Figz\fNIRS_Biv_HbR_ColorMap.png"

# Create versions with and without border
# create_legend(fnirs_hbr_colors, "Subnetworks", fnirs_hbr_output_base, with_border=True)
create_legend(fnirs_hbr_colors, "Subnetworks", fnirs_hbr_output_base, with_border=False)


# for partial


fmri_par_colors = [
    ([0.55, 0.0, 0.0], r'1$^{st}$: DMN'),
    ([0.10, 0.70, 0.20], r'2$^{nd}$: CEN'),
    ([0.20, 0.40, 0.90], r'3$^{rd}$: DMN'),
    ([0.98, 0.85, 0.20], r'4$^{th}$: d-AN'),
    ([0.80, 0.20, 0.80], r'5$^{th}$: SMN'),
    ([0.10, 0.80, 0.80], r'6$^{th}$: d-AN, p-VN, c-VN'),
    ([1.00, 0.50, 0.00], r'7$^{th}$: p-VN')
]

fmri_par_output_base = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\99_figures_Gimp\BrainNet Figz\fMRI_Par_ColorMap.png"
# Create versions with and without border
# create_legend(fmri_par_colors, "fMRI Subnetworks", fmri_par_output_base, with_border=True)
create_legend(fmri_par_colors, "fMRI Subnetworks", fmri_par_output_base, with_border=False)

fnirs_hbo_par_colors = [
    ([0.55, 0.0, 0.0], r'1$^{st}$'), # vivid red
    ([0.10, 0.70, 0.20], r'2$^{nd}$'), # emerald green
    ([0.80, 0.20, 0.80], r'3$^{rd}$'), # violet
    ([1.00, 0.50, 0.00], r'5$^{th}$'), # orange
    # ([0.20, 0.40, 0.90], r'3$^{rd}$'), # royal blue
    ([0.98, 0.85, 0.20], r'6$^{th}$'), # gold yellow

    # ([0.10, 0.80, 0.80], r'6$^{th}$'), # turquoise

]

fnirs_hbo_par_output_base = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\99_figures_Gimp\BrainNet Figz\fNIRS_Par_HbO_ColorMap.png"
# Create versions with and without border
# create_legend(fnirs_hbo_par_colors, "Subnetworks", fnirs_hbo
create_legend(fnirs_hbo_par_colors, "Subnetworks", fnirs_hbo_par_output_base, with_border=False)

fnirs_hbr_par_colors = [
    ([0.55, 0.0, 0.0], r'1$^{st}$'), # vivid red
    ([0.10, 0.70, 0.20], r'2$^{nd}$'), # emerald green
    ([0.80, 0.20, 0.80], r'4$^{th}$'), # violet
    ([1.00, 0.50, 0.00], r'7$^{th}$'), # orange
    # ([0.20, 0.40, 0.90], r'3$^{rd}$'), # royal blue
    # ([0.98, 0.85, 0.20], r'$^{th}$'), # gold yellow

    # ([0.10, 0.80, 0.80], r'6$^{th}$'), # turquoise
]
fnirs_hbr_par_output_base = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\99_figures_Gimp\BrainNet Figz\fNIRS_Par_HbR_ColorMap.png"
# Create versions with and without border
# create_legend(fnirs_hbr_par_colors, "Subnetworks", fnirs_hbr_par_output_base, with_border=True)
create_legend(fnirs_hbr_par_colors, "Subnetworks", fnirs_hbr_par_output_base, with_border=False)
