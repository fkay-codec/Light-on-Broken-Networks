"""
To visualize the TOST results across equivalence bounds, 
we plotted the percentage of significant edges as a function of Fisher-z transformed 
bounds. To obtain a smooth representation, we interpolated the data using a monotone
cubic spline (PCHIP). The equivalence bound corresponding to 50% significant edges was 
estimated by inverting the smoothed curve using linear interpolation.

This procedure allows for an approximate determination of the equivalence bound at which 
half of the edges show significant equivalence, accounting for the discrete sampling of 
tested bounds and reducing the effect of noise in the raw percentages


For dummies like me:

You have a bunch of x-values (bounds) and y-values (percent significant edges). 
Your y-values only exist at certain x-values, like steps: 0.01, 0.02, 0.03…

Now, you want to know where y reaches 50%.

Without interpolation: 
    you just look at your steps and pick the first x where y is ≥50%. 
    This is easy but not very exact.

With interpolation: 
    you draw a straight line (or smooth curve) between the two points around 50% and 
    see exactly where it hits 50%. This is more precise and nicer 

Example: If at bound 0.10 you have 48% and at 0.11 you have 53%, the real “50% bound” is between 0.10 and 0.11, not exactly at either. Interpolation finds that exact spot.

equiv bounds
 [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1  0.11 0.12 0.13 0.14
 0.15 0.16 0.17 0.18 0.19 0.2  0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28
 0.29 0.3  0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 ]
equiv bounds (z)
 [0.01   0.02   0.03   0.04   0.05   0.0601 0.0701 0.0802 0.0902 0.1003
 0.1104 0.1206 0.1307 0.1409 0.1511 0.1614 0.1717 0.182  0.1923 0.2027
 0.2132 0.2237 0.2342 0.2448 0.2554 0.2661 0.2769 0.2877 0.2986 0.3095
 0.3205 0.3316 0.3428 0.3541 0.3654 0.3769 0.3884 0.4001 0.4118 0.4236]
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, interp1d
from scipy.interpolate import PchipInterpolator
import pandas as pd


def plot_tost_curve(equiv_bounds_z, y_axis, title="TOST Across Fisher-z Equivalence Bounds",save_path=None):
    """
    Plot TOST results across Fisher-z equivalence bounds.

    Parameters
    ----------
    equiv_bounds_z : array-like
        Equivalence bounds in Fisher-z units (x-axis)
    y_axis : array-like
        Percentage of significant edges for each bound (y-axis)
    title : str
        Plot title
    """
    # Ensure numpy arrays
    equiv_bounds_z = np.asarray(equiv_bounds_z)
    y_axis = np.asarray(y_axis)
    
    # Clip y values to 0-100%
    y_axis = np.clip(y_axis, 0, 100)

    # # --- Smooth curve ---
    # x_smooth = np.linspace(equiv_bounds_z.min(), equiv_bounds_z.max(), 400)
    # y_smooth = make_interp_spline(equiv_bounds_z, y_axis)(x_smooth)

    # Smooth the curve
    x_smooth = np.linspace(equiv_bounds_z.min(), equiv_bounds_z.max(), 400)
    y_smooth = PchipInterpolator(equiv_bounds_z, y_axis)(x_smooth)

    # Compute 50% vertical bound using smoothed curve
    f_interp = interp1d(y_smooth, x_smooth, bounds_error=False, fill_value="extrapolate")
    bound_50 = f_interp(50)

    # # --- Compute 50% vertical bound ---
    # f_interp = interp1d(y_axis, equiv_bounds_z, fill_value="extrapolate")
    # bound_50 = f_interp(50)  # approximate bound where 50% edges significant

    # --- Plot ---
    plt.figure(figsize=(6.5, 5))
    plt.plot(x_smooth, y_smooth, color='tab:blue', linewidth=2.2, label='TOST (FDR)')
    # plt.text(0.02, 0.02, '*FDR-corrected', transform=plt.gca().transAxes, 
    #         fontsize=8, verticalalignment='bottom', style='italic')    
    plt.scatter(equiv_bounds_z, y_axis, alpha=0.0, color='tab:blue', edgecolor='none', s=0)#, zorder=3)

    # 50% lines
    plt.axhline(50, color='gray', linestyle='--', linewidth=1.3)
    plt.text(equiv_bounds_z.max()*0.97, 52, '50% threshold', color='gray', fontsize=10, ha='right')
    plt.axvline(bound_50, color='maroon', linestyle='--', linewidth=1, alpha=0.7)
    # plt.text(bound_50+0.005, 10, f'50%: z≈{bound_50:.3f} ⇔ r≈{np.tanh(bound_50):.3f}', color='red', fontsize=7.5, rotation=90, va='bottom')

    # --- Labels ---
    plt.xlabel("Equivalence Bound (Fisher z)", fontsize=12)
    plt.ylabel("Percentage of Significant Edges (%)", fontsize=12)
    plt.title(title, fontsize=18, weight='bold')

    # --- Axis formatting ---
    plt.xlim(equiv_bounds_z.min(), equiv_bounds_z.max())
    plt.ylim(0, 105)
    # plt.xticks(equiv_bounds_z[::5], [f"{z:.3f}" for z in equiv_bounds_z[::5]], rotation=90, fontsize=4)  # every 4th tick
    # plt.xticks(equiv_bounds_z, [f"{z:.3f}" for z in equiv_bounds_z], rotation=90, fontsize=2)  # show all x-ticks, smaller font
    # All ticks
    plt.xticks(equiv_bounds_z, 
            [f"{z:.3f}" if i % 5 == 0 else '' for i, z in enumerate(equiv_bounds_z)],
            rotation=90, 
            fontsize=9)  # adjust font as needed
    #alternative 50$ annotation
    # Add an arrow pointing to the 50% bound
    plt.annotate(f'z={bound_50:.3f} \nr={np.tanh(bound_50):.3f}',
                xy=(bound_50, 50), 
                xytext=(bound_50 + 0.03, 35),
                arrowprops=dict(arrowstyle='->', color='maroon', lw=1.2),
                fontsize=12, color='maroon', weight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='maroon', alpha=0.8))

    plt.yticks(np.arange(0, 110, 10))
    # plt.tick_params(axis='both', which='major', labelsize=10) # comment this if you play around with ticks
    # ...existing code...
    #! Not needed anymore
    # # Add shaded region below 50% threshold
    # plt.axhspan(0, 50, alpha=0.15, color='lightcoral', label='< 50% Significant', zorder=0)
    # plt.axhspan(50, 100, alpha=0.15, color='lightgreen', label='≥ 50% Significant', zorder=0)

    # --- Grid and legend ---
    # plt.grid(alpha=0.25)
    # plt.legend(frameon=False, fontsize=10, loc='lower right')

    plt.legend(frameon=True, fontsize=12, loc='lower right', #alignment='center',
            fancybox=True, shadow=True, framealpha=0.95)
    plt.tight_layout()
    # --- Save/Show  ---
    if save_path:
        plt.savefig(save_path, dpi=300)#, bbox_inches='tight')
    # plt.show()


# replicate the bounds
min_bound = 0.04
max_bound = 0.4
# create the equivalence bounds array that will increase by 0.01
equivalence_bounds = np.arange(min_bound, max_bound, 0.002).round(3)
# print("equiv bounds\n", equivalence_bounds)

# convert the equivalence bounds from r to z using Fisher transformation
equivalence_bounds_z = np.arctanh(equivalence_bounds).round(5)
# print("equiv bounds (z)\n", equivalence_bounds_z)
# print(len(equivalence_bounds_z))
# quit()

par_hbo_res_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\3_CrossModality_ttost\Partial_HbO"
par_hbr_res_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\3_CrossModality_ttost\Partial_HbR"

biv_hbo_res_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\3_CrossModality_ttost\Bivariate_HbO"
biv_hbr_res_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\3_CrossModality_ttost\Bivariate_HbR"

# lets create a dataframe to do that
tost_summary_partial = {
    "Equivalence_Bound_z": [],
    "HbO_Significant_Percentage": [],
    "HbR_Significant_Percentage": []
}

tost_summary_bivariate = {
    "Equivalence_Bound_z": [],
    "HbO_Significant_Percentage": [],
    "HbR_Significant_Percentage": []
}



# for each file in the path load the excel and count how many are significant with a function
def count_sigificant_edges_folder_path(folder_path):
    sig_counts=[]
    bounds=[]
    for file in os.listdir(folder_path):
        print(f"Processing file: {file}")
        df_fdr = pd.read_excel(os.path.join(folder_path, file), index_col=0, sheet_name='fdr_cor_p_value')
        # count total edges with numbers
        upper_ind = np.triu_indices(len(df_fdr), k=1)
        upper_values = df_fdr.values[upper_ind]
        sig_edges = np.sum(upper_values < 0.05)
        total_edge = upper_values.shape[0]
        print(f"Significant edges: {sig_edges} / {total_edge}")
        percentage = (sig_edges / total_edge) * 100
        sig_counts.append(percentage)
        bound_str = "0."+file.split('_')[4].replace(".xlsx", "") 
        bound_float = float(bound_str)
        print(f"Extracted bound: {bound_float} from filename.")
        bounds.append(float(bound_str))
    return bounds, sig_counts


par_hbo_bounds, par_hbo_counts = count_sigificant_edges_folder_path(par_hbo_res_path)  # test the function
# do a quick test with the bounds if they are the same as they should
# test_list = equivalence_bounds_z.tolist()
# if par_hbo_bounds == test_list:
#     print("Partial HbO bounds match expected equivalence bounds.")
# quit()
par_hbr_bounds, par_hbr_counts = count_sigificant_edges_folder_path(par_hbr_res_path)

biv_hbo_bounds, biv_hbo_counts = count_sigificant_edges_folder_path(biv_hbo_res_path)  # test the function
biv_hbr_bounds, biv_hbr_counts = count_sigificant_edges_folder_path(biv_hbr_res_path)

# now we can create the dataframes
tost_summary_partial['Equivalence_Bound_z'] = equivalence_bounds_z.tolist()
tost_summary_partial['HbO_Significant_Percentage'] = par_hbo_counts
tost_summary_partial['HbR_Significant_Percentage'] = par_hbr_counts

tost_summary_bivariate['Equivalence_Bound_z'] = equivalence_bounds_z.tolist()
tost_summary_bivariate['HbO_Significant_Percentage'] = biv_hbo_counts
tost_summary_bivariate['HbR_Significant_Percentage'] = biv_hbr_counts

# convert to dataframe
tost_summary_partial_df = pd.DataFrame(tost_summary_partial)
tost_summary_bivariate_df = pd.DataFrame(tost_summary_bivariate)

out_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\3_CrossModality_ttost\FiguresAltTicking"
if not os.path.exists(out_path):
    os.makedirs(out_path)
# plot the results
plot_tost_curve(
    tost_summary_partial_df["Equivalence_Bound_z"],
    tost_summary_partial_df["HbO_Significant_Percentage"],
    title="Edgewise Similarity of fMRI & fNIRS (HbO)\n(Partial Correlations)", #         "Partial Correlation: TOST Across Fisher-z Equivalence Bounds (fNIRS HbO vs fMRI)",
    save_path=os.path.join(out_path, "Partial_Correlation_TOST_HbO_vs_fMRI.png")
)

plot_tost_curve(
    tost_summary_partial_df["Equivalence_Bound_z"],
    tost_summary_partial_df["HbR_Significant_Percentage"],
    title="Edgewise Similarity of fMRI & fNIRS (HbR)\n(Partial Correlations)", #"Partial Correlation: TOST Across Fisher-z Equivalence Bounds (fNIRS HbR vs fMRI)",
    save_path=os.path.join(out_path, "Partial_Correlation_TOST_HbR_vs_fMRI.png")
)

plot_tost_curve(
    tost_summary_bivariate_df["Equivalence_Bound_z"],
    tost_summary_bivariate_df["HbO_Significant_Percentage"],
    title="Edgewise Similarity of fMRI & fNIRS (HbO)\n(Bivariate Correlations)", #"Bivariate Correlation: TOST Across Fisher-z Equivalence Bounds (fNIRS HbO vs fMRI)",
    save_path=os.path.join(out_path, "Bivariate_Correlation_TOST_HbO_vs_fMRI.png")
)

plot_tost_curve(
    tost_summary_bivariate_df["Equivalence_Bound_z"],
    tost_summary_bivariate_df["HbR_Significant_Percentage"],
    title="Edgewise Similarity of fMRI & fNIRS (HbR)\n(Bivariate Correlations)", #"Bivariate Correlation: TOST Across Fisher-z Equivalence Bounds (fNIRS HbR vs fMRI)",
    save_path=os.path.join(out_path, "Bivariate_Correlation_TOST_HbR_vs_fMRI.png")
)
