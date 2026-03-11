import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

def process_nodal_strength(nodal_df):
    """
    Process the dataframe

    Extract Nodal Str and then normalize it.

    At the end return a series with the normalized values
    
    """
    nodal_df = nodal_df.copy()
    nodal_df = nodal_df[["Positive Nodal Strength", "Negative Nodal Strength"]]

    
    # normalize based on Rubinov + Sporns 2011: [average positive] - average negative*(weighted by the proportion of negative to total)
    nodal_df["Normalized Str"] = (nodal_df["Positive Nodal Strength"]/81) - ((nodal_df["Negative Nodal Strength"]/81)*(nodal_df["Negative Nodal Strength"]/(nodal_df["Positive Nodal Strength"]+nodal_df["Negative Nodal Strength"])))

    return nodal_df["Normalized Str"]



# Custom function to load data depending on modality
def load_fmri_data_from_excel(file_path):
    df = pd.read_excel(file_path, sheet_name='GraphMetrics', index_col=0)
    return df

# Custom function to load data depending on modality
def load_fnirs_data_from_excel(file_path):
    hbo_df = pd.read_excel(file_path, sheet_name='HbO_GraphMetrics', index_col=0)
    hbr_df = pd.read_excel(file_path, sheet_name='HbR_GraphMetrics', index_col=0)
    return hbo_df, hbr_df

def compare_graph_metrics(fmri_df, fnirs_df, results_folder, name=""):
    """
    Correlate the graph metrics between fMRI and fNIRS dataframes.
    Save the results in a text file.
    Save the scatter plots as well.
    """

    fmri_df = fmri_df.copy()
    fnirs_df = fnirs_df.copy()
    # rename the index if they match
    fmri_df, fnirs_df = check_rename(fmri_df, fnirs_df)


    # clean irrelevant columns
    print(fmri_df.columns)
    columns_to_drop = [
        'Negative Assortativity',
        'Number of Positives','Number of Negatives',
    ]
    
    fmri_df=fmri_df.drop(columns=columns_to_drop, errors='raise')
    fnirs_df=fnirs_df.drop(columns=columns_to_drop, errors='raise')
    print(fmri_df.columns)
    print(fmri_df.head())


    results_df = pd.DataFrame(columns=[
        'Graph Metric', 'Pearson r', 'Pearson p', 'Spearman rho', 'Spearman p',
        'Kendall tau', 'Kendall p', 'FMRI Shapiro p', 'fNIRS Shapiro p',
        'Status', 'Notes'
    ])
    columns = fmri_df.columns
    for column in columns:
        if fmri_df[column].name != fnirs_df[column].name:
            print(f"Column names do not match: {fmri_df[column].name} vs {fnirs_df[column].name}... Exiting...")
            quit()

        fmri_values = fmri_df[column].values
        fnirs_values = fnirs_df[column].values
        print(f"Comparing Graph Metric: {column}...")
        # run stats while capturing warnings/exceptions and record status
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                _, shapiro_fmri_p = stats.shapiro(fmri_values)
                _, shapiro_fnirs_p = stats.shapiro(fnirs_values)
                r, p = stats.pearsonr(fmri_values, fnirs_values)
                rho, rho_p = stats.spearmanr(fmri_values, fnirs_values)
                tau, tau_p = stats.kendalltau(fmri_values, fnirs_values)
                warn_msgs = "; ".join(str(x.message) for x in w) if w else ""
            status = "ok" if not warn_msgs else "warning"
            notes = warn_msgs
        except Exception as e:
            # on error, record NaNs and the error message
            shapiro_fmri_p = np.nan
            shapiro_fnirs_p = np.nan
            r = p = rho = rho_p = tau = tau_p = np.nan
            status = "error"
            notes = str(e)

        # store the results in a row series and add it to the results df
        row_series = pd.Series({
            'Graph Metric': column,
            'Pearson r': r, 'Pearson p': p,
            'Spearman rho': rho, 'Spearman p': rho_p,
            'Kendall tau': tau, 'Kendall p': tau_p,
            'FMRI Shapiro p': shapiro_fmri_p,
            'fNIRS Shapiro p': shapiro_fnirs_p,
            'Status': status,
            'Notes': notes
        })
        results_df = pd.concat([results_df, row_series.to_frame().T], ignore_index=True)


    results_name = f"GraphMetricCorrelation_between_fMRI_and_fNIRS_{name}.xlsx"
    with pd.ExcelWriter(os.path.join(results_folder, results_name)) as writer:
        results_df.to_excel(writer, index=False, sheet_name='GraphMetric_Correlation')
    print(f"Results saved to {os.path.join(results_folder, results_name)}")

    # print(name)

    return 0

def check_rename(fmri_df, fnirs_df):
    """
    check if the column names match between the two dataframes: ROI1_S1D1 and S1D1
    if they match make them the same by just replacing with the fMRI one
    """
    fmri_df = fmri_df.copy()
    fnirs_df = fnirs_df.copy()

    index_if_match = fmri_df.index.tolist()
    # replace '_' with '/'
    index_if_match = [s.replace("_", "/") for s in index_if_match]
    
    # now in the fmri_df remove the prefix and then check if the columns match
    fmri_df.index = [s.split("_")[-1] for s in fmri_df.index.tolist()]
    # remove the "_" from the fnirs index as well
    fnirs_df.index = [s.replace("_","") for s in fnirs_df.index.tolist()]

    # now itterate through the fmri index and check if they match with the fnirs index
    if fmri_df.index.tolist() == fnirs_df.index.tolist():
        print("Index names match. Renaming fmri index and fnirs index.")
    else:
        print("Index names do not match... Exiting... ")
        quit()

    # rename both to the index_if_match
    fmri_df.index = index_if_match
    fnirs_df.index = index_if_match
    return fmri_df, fnirs_df

def plot_correlation_scatter_sns(fmri_values=None, fnirs_values=None, statistic=None, p_value=None, title=None, hemetype=None, correlation_type=None):
    """
    Plot the correlation between fMRI and fNIRS values using Seaborn's lmplot.
    Includes regression line and 95% confidence interval.
    """
    if correlation_type is None:
        raise ValueError("At least one correlation coefficient must be provided & specify its type.")

    if correlation_type=="Pearson":
        statistic_name = "r"
    elif correlation_type=="Spearman":
        statistic_name = "rho"
    elif correlation_type=="Kendall":
        statistic_name = "tau"

    # Create a tidy DataFrame for Seaborn
    df = pd.DataFrame({
        'fMRI values': fmri_values,
        'fNIRS values': fnirs_values,
        })
    
    # Create the lmplot (scatter + regression line)
    g = sns.lmplot(
        data=df,
        x='fMRI values',
        y='fNIRS values',
        ci=95,          # show 95% confidence interval
        scatter_kws={'alpha': 0.5, 'color': 'black', 's': 50},
        line_kws={'color': 'slateblue', 'lw': 2},
        height=8,
        aspect=1
    )
    
    # Customize the plot
    ax = g.ax
    ax.set_xlabel("fMRI values", fontsize=14)
    ax.set_ylabel("fNIRS values", fontsize=14)
    ax.set_title(
        f'{title}\n{correlation_type} Correlation between fMRI and fNIRS ({hemetype})',
        fontsize=16, fontweight='bold'
    )
    ax.grid(True, alpha=0.3)
    
    # Create legend entries
    # Bold r and p if p < 0.05
    if 0.01 <= p_value < 0.05:
        stats_label = f"{statistic_name} = {statistic:.2f}, p = {p_value:.2e}*"
    elif 0.001<= p_value < 0.01:
        stats_label = f"{statistic_name} = {statistic:.2f}, p = {p_value:.2e}**"
    elif p_value < 0.001:
        stats_label = f"{statistic_name} = {statistic:.2f}, p = {p_value:.2e}***"
    else:
        stats_label = f"{statistic_name} = {statistic:.2f}, p = {p_value:.2e}"
    line_label = "Fit line (95% CI)"
    # Create handles: one invisible (for text), one line for the fit
    stats_handle = plt.Line2D([], [], color='none', label=stats_label)  # invisible handle
    line_handle = plt.Line2D([], [], color='grey', linewidth=2, label=line_label)

    # Add both to the legend
    ax.legend(handles=[stats_handle, line_handle], fontsize=10, loc='lower right', frameon=True)
    plt.tight_layout()
    return g

def get_simplified_title(column_name):
    lower = column_name.lower()
    
    # Check for Nodal Path Length
    if "nodal path" in lower:
        if "3" in lower:
            return "Nodal Path Length 3"
        elif "4" in lower:
            return "Nodal Path Length 4"
    
    # Check for Nodal Betweenness
    elif "nodal betweenness" in lower:
        if "3" in lower:
            return "Nodal Betweenness 3"
        elif "4" in lower:
            return "Nodal Betweenness 4"
    
    # Otherwise, just return the original name
    return column_name

# set the output paths for results
results_master = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\2_CrossModalityGroupAnalysis\2_GraphMetricsComparison"
if not os.path.exists(results_master):
    os.makedirs(results_master)

# load the folders with the group graph metrics
fmri_folder_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\7_GraphMetricsGroupLevel_Updated"
fnirs_folder_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\5_GraphMetricsGroupLevel_Updated"

# extract the excel files from the folders
fmri_files = [f for f in os.listdir(fmri_folder_path) if f.endswith('.xlsx')]
fnirs_files = [f for f in os.listdir(fnirs_folder_path) if f.endswith('.xlsx')]

# print for verification
print("fMRI files:", fmri_files)
print("fNIRS files:", fnirs_files)

for fmri_file in fmri_files:
    fmri_df = pd.DataFrame()
    fnirs_df_hbo = pd.DataFrame()
    fnirs_df_hbr = pd.DataFrame()

    fmri_df = load_fmri_data_from_excel(os.path.join(fmri_folder_path, fmri_file))
    # add nodal strength normalized
    # fmri_df["Normalized Nodal Strength"] = process_nodal_strength(fmri_df)


    print(f"Loaded fMRI data from: \n\t{fmri_file}")
    # take the equivelent fnirs file
    i=0
    while fnirs_files[i].replace("fNIRS_","") != fmri_file.replace("fMRI_",""):
        print("Current fnirs file:", fnirs_files[i])
        i += 1
    print(f"Matched fnirs file: \n\t{fnirs_files[i]}")
    print(f"Computing comparisons...\n\n")

    # load the fnirs data (HbO and HbR)
    fnirs_df_hbo, fnirs_df_hbr = load_fnirs_data_from_excel(os.path.join(fnirs_folder_path, fnirs_files[i]))
    # fnirs_df_hbo["Normalized Nodal Strength"] = process_nodal_strength(fnirs_df_hbo)
    # fnirs_df_hbr["Normalized Nodal Strength"] = process_nodal_strength(fnirs_df_hbr)
    # extract the name for saving
    name = fmri_file.split('_')[3]  # "partial" or "standard"
    if not ("partial" in name.lower() or "standard" in name.lower()):
        raise ValueError(f"Filename -{name}- does not contain 'partial' or 'standard' to extract comparison type... Error will occur when writting")
    print("Comparison type:", name)
    results_folder = os.path.join(results_master, name)
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    compare_graph_metrics(fmri_df, fnirs_df_hbo, results_folder, name=name+"_HbO")
    compare_graph_metrics(fmri_df, fnirs_df_hbr, results_folder, name=name+"_HbR")

