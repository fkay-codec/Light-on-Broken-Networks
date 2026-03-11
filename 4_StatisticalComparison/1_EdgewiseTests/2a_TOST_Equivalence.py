"""
TODO: invesre z transf or z transform the bounds before tost?
"""
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
import progressbar
from numba import njit
from statsmodels.stats.weightstats import ttost_ind 
import pingouin as pg
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, interp1d

def fmri_load_all_subjects(file_path_list):
    """Load all fMRI subjects' data from the given file paths.
    We store them in a dictionary."""
    fmri_all_subjs = {}
    for file in file_path_list:     
        data = load_fmri_matrix(file)
        subj_id = os.path.basename(file).split('_')[0]
        fmri_all_subjs[subj_id] = data
    return fmri_all_subjs

def fnirs_load_all_subjects(file_path_list):
    """Load all fNIRS subjects' data from the given file paths.
    We first load the HbO and HbR matrices, clean the column and index names to match the ones in fmri (basically remove an underscore '_'), and store them in a dictionary.

    Dictionary format:
    {'subj_id': {'HbO': dataframe, 'HbR': dataframe}}

    We can navigate the dictionary as follows:
    fnirs[subj_id]['HbO'] to get the HbO dataframe for a specific subject
    etc."""
    fnirs_all_subjs = {}
    for file in file_path_list:     
        hbo, hbr = load_fnirs_matrix(file)
        subj_id = os.path.basename(file).split('_')[0]
        fnirs_all_subjs[subj_id] = {'HbO': hbo, 'HbR': hbr}
    return fnirs_all_subjs

def load_fmri_matrix(data_path):
    """
    Input: data_path: str
        Path to the Excel file containing fMRI data.
    Returns the dataframe from the given Excel file path.
    """
    df = pd.read_excel(data_path, index_col=0)
    return df

def load_fnirs_matrix(data_path):
    """
    Input: data_path: str
        Path to the Excel file containing HbO and HbR data.

    Returns the HbO and HbR dataframes from the given Excel file path.
    """

    hbo = pd.read_excel(data_path, sheet_name='HbO', index_col=0)
    # and clean the column and index names
    hbo.columns = [col.replace("_", "") for col in hbo.columns]
    hbo.index = [idx.replace("_", "") for idx in hbo.index]

    hbr = pd.read_excel(data_path, sheet_name='HbR', index_col=0)
    hbr.columns = [col.replace("_", "") for col in hbr.columns]
    hbr.index = [idx.replace("_", "") for idx in hbr.index]
    return hbo, hbr

def column_order_ismatch(fnirs_dict, fmri_dict):
    """
    Check if all dataframes in fnirs_dict and fmri_dict have the same column/index order.
    fnirs_dict format: {'subj_id': {'HbO': df, 'HbR': df}}
    fmri_dict format:  {'subj_id': df}

    If all dataframes have the same column/index order no error
    If there is a mismatch then it exits the program...
    """
    reference_columns = None
    
    # derive reference_columns from fnirs (use HbO from first subject)
    for subj_id, modality_dict in fnirs_dict.items():
        # print(f"Checking fnirs subject: {subj_id}")
        # for new structure, modality_dict is always {'HbO': df, 'HbR': df}
        if not isinstance(modality_dict, dict):
            print(f"Error: fnirs_dict[{subj_id}] is not a dict of modalities")
            print(f"Exitting...")
            quit()
        
        # get HbO as reference (or first available modality)
        df = modality_dict.get('HbO')
        if df is None:
            df = modality_dict.get('HbR')
        if df is None:
            df = next(iter(modality_dict.values()))

        cols = df.columns.tolist()
        rows = df.index.tolist()
        
        # check symmetry
        if cols != rows:
            print(f"Warning: fnirs df {subj_id} has non-symmetric labels (cols!=rows).")
            print(f"Exitting...")
            quit()

        # establish reference from first subject
        if reference_columns is None:
            reference_columns = cols
            # print(f"Reference columns established from subject {subj_id}: {len(reference_columns)} regions")
        elif cols != reference_columns:
            print(f"Column order mismatch in fnirs for subject: {subj_id}")
            print(f"Exitting...")
            quit()
        
        # check that both HbO and HbR have the same column order
        for mod_name, df_mod in modality_dict.items():
            # print(f"  Checking modality: {mod_name}")
            cols_mod = df_mod.columns.tolist()
            rows_mod = df_mod.index.tolist()
            
            if cols_mod != reference_columns:
                print(f"Column order mismatch in fnirs ({mod_name}) for subject: {subj_id}")
                quit()
            
            if rows_mod != reference_columns:
                print(f"Row order mismatch in fnirs ({mod_name}) for subject: {subj_id}")
                quit()

    if reference_columns is None:
        print("No fnirs data found to derive reference columns.")
        return False

    # check fmri dict - normalize labels by removing ROI prefix
    for subj_id, df in fmri_dict.items():
        # normalize the fmri column/row labels (remove network prefix like "Network_01_")
        norm_cols = [c.split("_", 1)[-1] if "_" in c else c for c in df.columns.tolist()]
        norm_rows = [r.split("_", 1)[-1] if "_" in r else r for r in df.index.tolist()]
        
        # check symmetry
        if norm_cols != norm_rows:
            print(f"Warning: fmri df {subj_id} has non-symmetric labels after normalization.")
            print(f"Exitting...")
            quit()
        # check against reference
        if norm_cols != reference_columns:
            print(f"Column order mismatch in fmri for subject: {subj_id}")
            print(f"Expected: {reference_columns[:5]}... ({len(reference_columns)} total)")
            print(f"Found:    {norm_cols[:5]}... ({len(norm_cols)} total)")
            print(f"Exitting...")
            quit()

    print(f"All dataframes have matching column/index orders.")
    return 1

def fdr_on_df(df):
    upper_ind = np.triu_indices(len(df.columns), k=1)
    p_values_upper = df.values[upper_ind]
    p_values_upper = np.array(p_values_upper, dtype=float)
    # apply fdr correction
    p_values_fdr = stats.false_discovery_control(p_values_upper, method='bh')
    # fill it back to the matrix
    fdr_corrected_df = df.copy()
    fdr_corrected_df.values[upper_ind] = p_values_fdr
    fdr_corrected_df.values[np.tril_indices(len(fdr_corrected_df), k=0)] = np.nan
    return fdr_corrected_df

def plot_tost_curve(equiv_bounds_z, y_axis, title="TOST Across Fisher-z Equivalence Bounds"):
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

    # --- Smooth curve ---
    x_smooth = np.linspace(equiv_bounds_z.min(), equiv_bounds_z.max(), 400)
    y_smooth = make_interp_spline(equiv_bounds_z, y_axis)(x_smooth)

    # --- Compute 50% vertical bound ---
    f_interp = interp1d(y_axis, equiv_bounds_z, fill_value="extrapolate")
    bound_50 = f_interp(50)  # approximate bound where 50% edges significant

    # --- Plot ---
    plt.figure(figsize=(7, 5))
    plt.plot(x_smooth, y_smooth, color='tab:blue', linewidth=2.5, label='TOST Results')
    plt.scatter(equiv_bounds_z, y_axis, color='tab:blue', edgecolor='white', s=50, zorder=3)

    # 50% lines
    plt.axhline(50, color='gray', linestyle='--', linewidth=1.3)
    plt.text(equiv_bounds_z.max()*0.97, 52, '50% threshold', color='gray', fontsize=10, ha='right')
    plt.axvline(bound_50, color='maroon', linestyle='--', linewidth=1.5, alpha=0.7)
    plt.text(bound_50+0.005, 10, f'50% at z≈{bound_50:.3f}', color='red', fontsize=10, rotation=90, va='bottom')

    # --- Labels ---
    plt.xlabel("Equivalence Bound (Fisher z)", fontsize=12)
    plt.ylabel("Percentage of Significant Edges (%)", fontsize=12)
    plt.title(title, fontsize=13, weight='bold')

    # --- Axis formatting ---
    plt.xlim(equiv_bounds_z.min()-0.01, equiv_bounds_z.max()+0.01)
    plt.ylim(0, 105)
    plt.xticks(equiv_bounds_z[::4], [f"{z:.3f}" for z in equiv_bounds_z[::4]])  # every 4th tick
    plt.yticks(np.arange(0, 110, 10))
    plt.tick_params(axis='both', which='major', labelsize=10)

    # --- Grid and legend ---
    plt.grid(alpha=0.25)
    plt.legend(frameon=False, fontsize=10, loc='lower right')

    plt.tight_layout()
    plt.show()

# load the paths
fmri_folder_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\3_CorrelationMatrices"
fnirs_folder_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\1_CorrelationMatrices"

# construct the paths for bivariate and partial correlation folders
# for StandardCorr: folder name=1_StandardCor || for PartialCorr: folder 2_PartialCor
fmri_bivariate_folder_path = os.path.join(fmri_folder_path, '1_StandardCor')
fnirs_bivariate_folder_path = os.path.join(fnirs_folder_path, '1_StandardCor')

fmri_partial_folder_path = os.path.join(fmri_folder_path, '2_PartialCor')
fnirs_partial_folder_path = os.path.join(fnirs_folder_path, '2_PartialCor')

# load the path of all files
fmri_bivariate_files = [os.path.join(fmri_bivariate_folder_path, f) for f in os.listdir(fmri_bivariate_folder_path) if f.endswith('.xlsx')]
fnirs_bivariate_files = [os.path.join(fnirs_bivariate_folder_path, f) for f in os.listdir(fnirs_bivariate_folder_path) if f.endswith('.xlsx')] 

fmri_partial_files = [os.path.join(fmri_partial_folder_path, f) for f in os.listdir(fmri_partial_folder_path) if f.endswith('.xlsx')]
fnirs_partial_files = [os.path.join(fnirs_partial_folder_path, f) for f in os.listdir(fnirs_partial_folder_path) if f.endswith('.xlsx')]

# load the data in a dictionary
# fMRI format: {'subj_id': dataframe}
# fNIRS format: {'subj_id': {{'HbO': dataframe, 'HbR': dataframe}}
fmri_bivariate_data = fmri_load_all_subjects(fmri_bivariate_files)
fnirs_bivariate_data = fnirs_load_all_subjects(fnirs_bivariate_files)
fmri_partial_data = fmri_load_all_subjects(fmri_partial_files)
fnirs_partial_data = fnirs_load_all_subjects(fnirs_partial_files)

# Sanity Check: are all columns and index names the same between fNIRS and fMRI dataframes?
# Perform the check
column_order_ismatch(fnirs_bivariate_data, fmri_bivariate_data)
column_order_ismatch(fnirs_partial_data, fmri_partial_data)

# we need to create a dictionary to store the results now...
# because we are performing the TOST for a range of equivalence bounds, the dictionary will have the following format:
# results_partial = {
#     "HbO Results": {
#         bound1: {
#             "p_value" : p_value_df1,                  # channel x channel dataframe of p-values || one p value per edge across all subjs
#             "fdr_cor_p_value" : fdr_cor_p_value_df1,  # channel x channel dataframe of fdr corrected p-values || one p value per edge across all subjs
#             "t-value" : t_value_df1,                  # channel x channel dataframe of t-values || one t value per edge across all subjs
#         }
#         bound2: {
#             "p_value" : p_value_df2,
#             "fdr_cor_p_value" : fdr_cor_p_value_df2,
#             "t-value" : t_value_df2,
#         }
#         ...
#     },
#     "HbR Results": {
#         ... same structure as HbO ...
#     }
# }

# lets create a template dataframe to store the results. The template will have the index/column of the channel names of fmri.

fmri_subjects = list(fmri_bivariate_data.keys())
fmri_first_subj = fmri_bivariate_data[fmri_subjects[0]]
reference_column = fmri_first_subj.columns
reference_column = [c.replace("_", "/") for c in reference_column]  # normalize the column names for better visualization if needed later... I dont think that we are going to visualize though as a matrix

print("reference_column\n", len(reference_column),"\n" , reference_column)

# create a template dataframe
template_df = pd.DataFrame(index=reference_column, columns=reference_column)

print("template_df\n", template_df)

# the equivalence bounds to test are like this because our data come from correlations... so they are in the range -1 to 1, thus we choose small values
min_bound = 0.04
max_bound = 0.4

# create the equivalence bounds array that will increase by 0.01
equivalence_bounds = np.arange(min_bound, max_bound, 0.002).round(3)
print("equiv bounds\n", equivalence_bounds)

# convert the equivalence bounds from r to z using Fisher transformation
equivalence_bounds_z = np.arctanh(equivalence_bounds).round(5)
print("equiv bounds (z)\n", equivalence_bounds_z)

# create the results dictionary for partial
results_partial = {
    "HbO Results": {},
    "HbR Results": {}
}

results_bivariate = {
    "HbO Results": {},
    "HbR Results": {}
}

# perform the TOST for both partial/bivariate hbo/hbr and store the results in their respective dictionaries
for bound in equivalence_bounds_z:
    print(f"Performing TOST for equivalence bound: {bound} (z units)")
    # for each bound create a new dataframe to store the p-values for each edge
    res_par_df_hbo_p = template_df.copy()
    res_par_df_hbr_p = template_df.copy()
    res_biv_df_hbo_p = template_df.copy()
    res_biv_df_hbr_p = template_df.copy()

    # create dataframes to store the fdr corrected p-values
    res_par_df_hbo_fdr_p = template_df.copy()
    res_par_df_hbr_fdr_p = template_df.copy()
    res_biv_df_hbo_fdr_p = template_df.copy()
    res_biv_df_hbr_fdr_p = template_df.copy()

    for i in range(len(template_df.index)): # rows
        for j in range(i+1, len(template_df.columns)): # columns
            # get the edge values across all subjects for partial correlation
            fmri_par_edge_values = []
            fnirs_par_edge_values_hbo = []
            fnirs_par_edge_values_hbr = []
            for subj_id in fmri_partial_data.keys():
                fmri_par_edge_values.append(fmri_partial_data[subj_id].iloc[i, j])
            for subj_id in fnirs_partial_data.keys():
                fnirs_par_edge_values_hbo.append(fnirs_partial_data[subj_id]['HbO'].iloc[i, j])
                fnirs_par_edge_values_hbr.append(fnirs_partial_data[subj_id]['HbR'].iloc[i, j])
            # convert to numpy arrays
            fmri_par_edge_values = np.array(fmri_par_edge_values)
            fnirs_par_edge_values_hbo = np.array(fnirs_par_edge_values_hbo)
            fnirs_par_edge_values_hbr = np.array(fnirs_par_edge_values_hbr)

            # perform TOST for HbO
            # penguin way
            # it does Welch correction to adjust for unequal variance...
            par_tost_hbo = pg.tost(fmri_par_edge_values, fnirs_par_edge_values_hbo, bound = bound, paired=False, correction=True)
            par_tost_hbo_p = par_tost_hbo['pval'].values

            # for HbR
            par_tost_hbr = pg.tost(fmri_par_edge_values, fnirs_par_edge_values_hbr, bound = bound, paired=False, correction=True)
            par_tost_hbr_p = par_tost_hbr['pval'].values

            # store them in their respective dataframe
            res_par_df_hbo_p.iloc[i, j] = par_tost_hbo_p
            res_par_df_hbr_p.iloc[i, j] = par_tost_hbr_p

            # We do the same for bivariate now
            # get the edge values across all subjects for bivariate correlation
            fmri_biv_edge_values = []
            fnirs_biv_edge_values_hbo = []
            fnirs_biv_edge_values_hbr = []
            for subj_id in fmri_bivariate_data.keys():
                fmri_biv_edge_values.append(fmri_bivariate_data[subj_id].iloc[i, j])
            for subj_id in fnirs_bivariate_data.keys():
                fnirs_biv_edge_values_hbo.append(fnirs_bivariate_data[subj_id]['HbO'].iloc[i, j])
                fnirs_biv_edge_values_hbr.append(fnirs_bivariate_data[subj_id]['HbR'].iloc[i, j])
            # convert to numpy arrays
            fmri_biv_edge_values = np.array(fmri_biv_edge_values)
            fnirs_biv_edge_values_hbo = np.array(fnirs_biv_edge_values_hbo)
            fnirs_biv_edge_values_hbr = np.array(fnirs_biv_edge_values_hbr)
            # perform TOST for HbO
            biv_tost_hbo = pg.tost(fmri_biv_edge_values, fnirs_biv_edge_values_hbo, bound = bound, paired=False, correction=True)
            biv_tost_hbo_p = biv_tost_hbo['pval'].values
            # for HbR
            biv_tost_hbr = pg.tost(fmri_biv_edge_values, fnirs_biv_edge_values_hbr, bound = bound, paired=False, correction=True)
            biv_tost_hbr_p = biv_tost_hbr['pval'].values
            # store them in their respective dataframe
            res_biv_df_hbo_p.iloc[i, j] = biv_tost_hbo_p
            res_biv_df_hbr_p.iloc[i, j] = biv_tost_hbr_p
    # compute the FDR correction DFs now for each modality and correlation type
    # Partial HbO/HbR
    res_par_df_hbr_fdr_p = fdr_on_df(res_par_df_hbr_p)  
    res_par_df_hbo_fdr_p = fdr_on_df(res_par_df_hbo_p)

    # Bivariate HbO/HbR
    res_biv_df_hbr_fdr_p = fdr_on_df(res_biv_df_hbr_p)  
    res_biv_df_hbo_fdr_p = fdr_on_df(res_biv_df_hbo_p)

    # store the dataframes in the results dictionary
    # Partial
    results_partial["HbO Results"][bound] = {
        "p_value" : res_par_df_hbo_p,
        "fdr_cor_p_value" : res_par_df_hbo_fdr_p,
    }
    results_partial["HbR Results"][bound] = {
        "p_value" : res_par_df_hbr_p,
        "fdr_cor_p_value" : res_par_df_hbr_fdr_p,
    }
    # Bivariate
    results_bivariate["HbO Results"][bound] = {
        "p_value" : res_biv_df_hbo_p,
        "fdr_cor_p_value" : res_biv_df_hbo_fdr_p,
    }
    results_bivariate["HbR Results"][bound] = {
        "p_value" : res_biv_df_hbr_p,
        "fdr_cor_p_value" : res_biv_df_hbr_fdr_p,
    }
  
    
# save the results to excel files now
output_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\3_CrossModality_ttost"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Partial
for modality in results_partial.keys():  # HbO Results, HbR Results
    modality_folder = os.path.join(output_folder, f"Partial_{modality.split()[0]}")  # Partial_HbO, Partial_HbR
    if not os.path.exists(modality_folder):
        os.makedirs(modality_folder)
    for bound, result_dict in results_partial[modality].items():
        # create an excel writer
        bound_clean = str(bound).replace('.', '_')
        output_file_path = os.path.join(modality_folder, f"ttost_results_bound_{bound_clean}.xlsx")
        with pd.ExcelWriter(output_file_path) as writer:
            # save p_value dataframe
            result_dict["p_value"].to_excel(writer, sheet_name='p_value')
            # save fdr corrected p_value dataframe
            result_dict["fdr_cor_p_value"].to_excel(writer, sheet_name='fdr_cor_p_value')

# Bivariate
for modality in results_bivariate.keys():  # HbO Results, HbR Results
    modality_folder = os.path.join(output_folder, f"Bivariate_{modality.split()[0]}")  # Bivariate_HbO, Bivariate_HbR
    if not os.path.exists(modality_folder):
        os.makedirs(modality_folder)
    for bound, result_dict in results_bivariate[modality].items():
        # create an excel writer
        bound_clean = str(bound).replace('.', '_')
        output_file_path = os.path.join(modality_folder, f"ttost_results_bound_{bound_clean}.xlsx")
        with pd.ExcelWriter(output_file_path) as writer:
            # save p_value dataframe
            result_dict["p_value"].to_excel(writer, sheet_name='p_value')
            # save fdr corrected p_value dataframe
            result_dict["fdr_cor_p_value"].to_excel(writer, sheet_name='fdr_cor_p_value')
print("All results saved successfully.")


# Not the end... #! maybe the end because it takes to long to do all these stuff...

# Now that we have the results saved we can do some basic vis.

# For example: we could count for each bound the precentage of significant edges and plot on x-axis the bounds (z) and on y-axis the percentage (%)
# one plot for HbO partial vs fMRI and one for HbR partial vs fMRI || and the same for bivariate 

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

for bound in equivalence_bounds_z:
    # Partial
    tost_summary_partial["Equivalence_Bound_z"].append(bound)
    partial_hbo_fdr_df = results_partial["HbO Results"][bound]["fdr_cor_p_value"]
    partial_hbr_fdr_df = results_partial["HbR Results"][bound]["fdr_cor_p_value"]
    # count significant edges
    total_edges = np.triu_indices(len(partial_hbo_fdr_df), k=1)[0].shape[0]
    print("total edges:", total_edges)
    
    significant_edges_hbo = np.nansum(partial_hbo_fdr_df.values[np.triu_indices(len(partial_hbo_fdr_df), k=1)] < 0.05)
    significant_edges_hbr = np.nansum(partial_hbr_fdr_df.values[np.triu_indices(len(partial_hbr_fdr_df), k=1)] < 0.05)
    percent_significant_hbo = (significant_edges_hbo / total_edges) * 100
    percent_significant_hbr = (significant_edges_hbr / total_edges) * 100
    tost_summary_partial["HbO_Significant_Percentage"].append(percent_significant_hbo)
    tost_summary_partial["HbR_Significant_Percentage"].append(percent_significant_hbr)
    # Bivariate
    tost_summary_bivariate["Equivalence_Bound_z"].append(bound)
    bivariate_hbo_fdr_df = results_bivariate["HbO Results"][bound]["fdr_cor_p_value"]
    bivariate_hbr_fdr_df = results_bivariate["HbR Results"][bound]["fdr_cor_p_value"]
    # count significant edges
    significant_edges_hbo_biv = np.nansum(bivariate_hbo_fdr_df.values[np.triu_indices(len(bivariate_hbo_fdr_df), k=1)] < 0.05)
    significant_edges_hbr_biv = np.nansum(bivariate_hbr_fdr_df.values[np.triu_indices(len(bivariate_hbr_fdr_df), k=1)] < 0.05)
    percent_significant_hbo_biv = (significant_edges_hbo_biv / total_edges) * 100
    percent_significant_hbr_biv = (significant_edges_hbr_biv / total_edges) * 100
    tost_summary_bivariate["HbO_Significant_Percentage"].append(percent_significant_hbo_biv)
    tost_summary_bivariate["HbR_Significant_Percentage"].append(percent_significant_hbr_biv)
# convert to dataframe
tost_summary_partial_df = pd.DataFrame(tost_summary_partial)
tost_summary_bivariate_df = pd.DataFrame(tost_summary_bivariate)

# plot the results
plot_tost_curve(
    tost_summary_partial_df["Equivalence_Bound_z"],
    tost_summary_partial_df["HbO_Significant_Percentage"],
    title="Partial Correlation: TOST Across Fisher-z Equivalence Bounds (fNIRS HbO vs fMRI)"
)

plot_tost_curve(
    tost_summary_partial_df["Equivalence_Bound_z"],
    tost_summary_partial_df["HbR_Significant_Percentage"],
    title="Partial Correlation: TOST Across Fisher-z Equivalence Bounds (fNIRS HbR vs fMRI)"
)

plot_tost_curve(
    tost_summary_bivariate_df["Equivalence_Bound_z"],
    tost_summary_bivariate_df["HbO_Significant_Percentage"],
    title="Bivariate Correlation: TOST Across Fisher-z Equivalence Bounds (fNIRS HbO vs fMRI)"
)

plot_tost_curve(
    tost_summary_bivariate_df["Equivalence_Bound_z"],
    tost_summary_bivariate_df["HbR_Significant_Percentage"],
    title="Bivariate Correlation: TOST Across Fisher-z Equivalence Bounds (fNIRS HbR vs fMRI)"
)
