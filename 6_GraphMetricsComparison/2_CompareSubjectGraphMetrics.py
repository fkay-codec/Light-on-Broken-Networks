import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

"""
we are going to perform an independent samples ttest on the nodal metrics.
for node across participants compare the fnirs vs fmri difference in the graph metric and save an excel
"""

txt = [] # this is a global list to print the summary at the end

# now for each combined graph metric we perform the ttest across subjects for each node
def perform_ttest_dictionary(fnirs_dict, fmri_dict):
    """
    Input:
        two dictionaries with the same keys representing graph metrics
    Return:
        a dictionary with the ttest results for each graph metric
    """
    def check_if_same_index(fmri, fnirs):
        fnirs = fnirs.copy()
        fmri = fmri.copy()
        # remove ROIXX from fmri
        fmri.index = [f.split('_')[-1] for f in fmri.index]

        if len(fmri_df.index) != len(fnirs_df.index):
            print("The number of nodes in fMRI and fNIRS do not match!")
            return False
        for node1, node2 in zip(fmri.index, fnirs.index):
            if node1 != node2:
                print(f"Node mismatch: fMRI node {node1} != fNIRS node {node2}")
                return False
        return True

    ttest_results_dict = {}
    i = 0
    for metric in fmri_dict.keys():     
        print(f"    Performing ttest for metric: {metric}")
        fmri_df = fmri_dict[metric].copy()
        fnirs_df = fnirs_dict[metric].copy()
        # print(fnirs_df.head())
        # quit()
        # check if the rows (nodes) are the same
        if i == 0:
            if check_if_same_index(fmri=fmri_df, fnirs=fnirs_df):
                print("Indices match, proceeding with ttest...")
            else:
                print("Indices do not match, cannot perform ttest.")
                quit()
        fnirs_df.index = fmri_df.index  # Ensure indices match for ttest
        
        # print(fmri_df.head())
        # print(fnirs_df.head())
        # quit()
        for node in fmri_df.index:
            # print(fnirs_df.loc[node])
            # quit()
            # print(f"Performing t-test for Node: {node}")
            fmri_values = fmri_df.loc[node].values.astype(float)
            # print(fmri_values)

            fnirs_values = fnirs_df.loc[node].values.astype(float)
            # print(fnirs_values)
            std1 = stats.tstd(fmri_values)
            # print(f"std1: {std1}")
            # quit()
            std2 = stats.tstd(fnirs_values)
            # print(f"std2: {std2}")
            # quit()
            std_pooled = np.sqrt(((std1 **2) + (std2 **2)) /2)
            fmri_mean = np.mean(fmri_values)
            fnirs_mean = np.mean(fnirs_values)
            # safety: avoid division by zero when pooled std is zero
            if std_pooled == 0 or not np.isfinite(std_pooled):
                cohens_d = np.nan
                note = 'std_pooled_zero_or_nan'
            else:
                cohens_d =  (fmri_mean - fnirs_mean) / std_pooled
                note = ''
            t_stat, p_value = stats.ttest_ind(fmri_values, fnirs_values, alternative='two-sided')
            # print(f"Cohen's d for Node {node}, Metric {metric}: {cohens_d}")
            # print(f"mean fMRI: {fmri_mean}, mean fNIRS: {fnirs_mean}")
            # print(f"t-statistic: {t_stat}, p-value: {p_value}")
            # # quit()

            if metric not in ttest_results_dict:
                ttest_results_dict[metric] = pd.DataFrame(columns=['t-stat', 'p-value', 'fdr-p', 'fmri_mean', 'fnirs_mean', 'std_pooled', 'cohens_d', 'Notes'], index = fmri_df.index)
            # print(ttest_results_dict[metric].head())
            # quit()
            ttest_results_dict[metric].loc[node] = {
                    't-stat': t_stat,
                    'p-value': p_value,
                    'fdr-p': np.nan,  # placeholder, will be filled after FDR correction
                    'fmri_mean': fmri_mean,
                    'fnirs_mean': fnirs_mean,
                    'std_pooled': std_pooled,
                    'cohens_d': cohens_d,
                    'Notes': note
                }
            # print(f"Node: {node}, t-stat: {t_stat}, p-value: {p_value}")
        # print(ttest_results_dict[metric].head())
        # apply fdr correction
        # sanitize p-values before FDR: keep only finite values in [0,1]
        pvals = ttest_results_dict[metric]['p-value'].astype(float).to_numpy()
        valid_mask = np.isfinite(pvals) & (pvals >= 0.0) & (pvals <= 1.0)
        fdr_array = np.full_like(pvals, np.nan, dtype=float)
        if valid_mask.any():
            try:
                corrected = stats.false_discovery_control(pvals[valid_mask], method='bh')
                fdr_array[valid_mask] = corrected
            except Exception:
                # if FDR fails for any reason, leave fdr as NaN
                fdr_array[valid_mask] = np.nan
        ttest_results_dict[metric]['fdr-p'] = fdr_array
        # print(ttest_results_dict[metric].head())

        # quit()
        i += 1


        #! silly stuff that are against compute science but whatever
        global txt  # Declare that we use the global variable
        count_significant = (ttest_results_dict[metric]['fdr-p'] < 0.05).sum()
        # print(f"    Number of significant nodes after FDR correction for metric {metric}: {count_significant}")
        txt.append(f"       {metric}: {round((count_significant/82)*100, 2)}%")  # Store result in global txt
    return ttest_results_dict


# Custom function to load data depending on modality
def load_fmri_data_from_excel(file_path):
    df = pd.read_excel(file_path, sheet_name='GraphMetrics', index_col=0)
    return df
# Custom function to load data depending on modality
def load_fnirs_data_from_excel(file_path):
    hbo_df = pd.read_excel(file_path, sheet_name='HbO_GraphMetrics', index_col=0)
    hbr_df = pd.read_excel(file_path, sheet_name='HbR_GraphMetrics', index_col=0)
    return hbo_df, hbr_df

# we could combine all the graph metrics into a single dataframe for easier handling... it will increase computational time however but it will be cleaner
def clean_from_irrelevant_graph_metrics(df):
    df = df.copy()
    # Remove columns that are not graph metrics

    relevant_metrics = ['Net Nodal Strength', 'Normalized Nodal Strength','Positive Assortativity', 'Local Efficiency','Nodal Path Length [1-|z|]', 'Betweenness Centrality [1-|z|]']
    cleaned_df = df[relevant_metrics]
    # rename the column names for easier identification

    new_names = ['Net NodStr', 'Norm NodStr', 'Pos Assort', 'Loc Eff', 'Nod PathLen', 'Nod Betweenness']
    cleaned_df.columns = new_names
    return cleaned_df

def add_subj_name_to_columns(df, subj_name):
    df = df.copy()
    df.columns = [f"{subj_name}_{col}" for col in df.columns]
    return df

# perform it for bivariate first
def combine_fnirs_data_to_dict(fnirs_folder_path):
    
    """
    We are reading the excel files from the fnirs folder and then we create a dataframe for HbO and HbR with each column representing a subject and each row the nodes
    then we place it in a dictionary.
    the dictionary will look like this:
    dictionary = {
        'GraphMetric1': Combined_hbo_data_for_that_metric,
        'GraphMetric2': Combined_hbr_data_for_that_metric,
        ...
    }
    We will return two dictionaries, one for HbO and one for HbR   
    """
    fnirs_files = [f for f in os.listdir(fnirs_folder_path) if f.endswith('.xlsx')]
    combined_hbo_data = {
        'Net Nodal Strength': pd.DataFrame(),
        'Norm Nodal Strength': pd.DataFrame(),
        'Positive Assortativity': pd.DataFrame(),
        'Local Efficiency': pd.DataFrame(),
        'Nodal Path Length': pd.DataFrame(),
        'Nodal Betweenness': pd.DataFrame()
    }
    combined_hbr_data = {
        'Net Nodal Strength': pd.DataFrame(),
        'Norm Nodal Strength': pd.DataFrame(),
        'Positive Assortativity': pd.DataFrame(),
        'Local Efficiency': pd.DataFrame(),
        'Nodal Path Length': pd.DataFrame(),
        'Nodal Betweenness': pd.DataFrame()
    }

    for fnirs_subj in fnirs_files:

        fnirs_subj_path = os.path.join(fnirs_folder_path, fnirs_subj)
        hbo_df, hbr_df = load_fnirs_data_from_excel(fnirs_subj_path)

        clean_hbo_df = clean_from_irrelevant_graph_metrics(hbo_df)
        clean_hbr_df = clean_from_irrelevant_graph_metrics(hbr_df)
        # remove the '_' from the index names
        clean_hbo_df.index = clean_hbo_df.index.str.replace('_', '')
        clean_hbr_df.index = clean_hbr_df.index.str.replace('_', '')


        # extract name of subj
        subj_name = fnirs_subj.split('_')[0]

        combined_hbo_data['Net Nodal Strength'] = pd.concat([combined_hbo_data['Net Nodal Strength'], clean_hbo_df['Net NodStr'].rename(subj_name + "_Net NodStr")], axis=1)
        combined_hbo_data['Norm Nodal Strength'] = pd.concat([combined_hbo_data['Norm Nodal Strength'], clean_hbo_df['Norm NodStr'].rename(subj_name + "_Norm NodStr")], axis=1)
        combined_hbo_data['Positive Assortativity'] = pd.concat([combined_hbo_data['Positive Assortativity'], clean_hbo_df['Pos Assort'].rename(subj_name + "_Pos Assort")], axis=1)
        combined_hbo_data['Local Efficiency'] = pd.concat([combined_hbo_data['Local Efficiency'], clean_hbo_df['Loc Eff'].rename(subj_name + "_Loc Eff")], axis=1)
        combined_hbo_data['Nodal Path Length'] = pd.concat([combined_hbo_data['Nodal Path Length'], clean_hbo_df['Nod PathLen'].rename(subj_name + "_Nod PathLen")], axis=1)
        combined_hbo_data['Nodal Betweenness'] = pd.concat([combined_hbo_data['Nodal Betweenness'], clean_hbo_df['Nod Betweenness'].rename(subj_name + "_Nod Betweenness")], axis=1)    

        combined_hbr_data['Net Nodal Strength'] = pd.concat([combined_hbr_data['Net Nodal Strength'], clean_hbr_df['Net NodStr'].rename(subj_name + "_Net NodStr")], axis=1)
        combined_hbr_data['Norm Nodal Strength'] = pd.concat([combined_hbr_data['Norm Nodal Strength'], clean_hbr_df['Norm NodStr'].rename(subj_name + "_Norm NodStr")], axis=1)
        combined_hbr_data['Positive Assortativity'] = pd.concat([combined_hbr_data['Positive Assortativity'], clean_hbr_df['Pos Assort'].rename(subj_name + "_Pos Assort")], axis=1)
        combined_hbr_data['Local Efficiency'] = pd.concat([combined_hbr_data['Local Efficiency'], clean_hbr_df['Loc Eff'].rename(subj_name + "_Loc Eff")], axis=1)
        combined_hbr_data['Nodal Path Length'] = pd.concat([combined_hbr_data['Nodal Path Length'], clean_hbr_df['Nod PathLen'].rename(subj_name + "_Nod PathLen")], axis=1)
        combined_hbr_data['Nodal Betweenness'] = pd.concat([combined_hbr_data['Nodal Betweenness'], clean_hbr_df['Nod Betweenness'].rename(subj_name + "_Nod Betweenness")], axis=1)  

    return combined_hbo_data, combined_hbr_data

def combine_fmri_data_to_dict(fmri_folder_path):
    
    """
    We are reading the excel files from the fmri folder and then we create a dataframe for fmri with each column representing a subject and each row the nodes
    then we place it in a dictionary.
    the dictionary will look like this:
    dictionary = {
        'GraphMetric1': Combined_fmri_data_for_that_metric,
        'GraphMetric2': Combined_fmri_data_for_that_metric,
        ...
    }
    We will return the dictionary  
    """
    fmri_files = [f for f in os.listdir(fmri_folder_path) if f.endswith('.xlsx')]
    combined_fmri_data = {
        'Net Nodal Strength': pd.DataFrame(),
        'Norm Nodal Strength': pd.DataFrame(),
        'Positive Assortativity': pd.DataFrame(),
        'Local Efficiency': pd.DataFrame(),
        'Nodal Path Length': pd.DataFrame(),
        'Nodal Betweenness': pd.DataFrame(),

    }

    for fmri_subj in fmri_files:

        fmri_subj_path = os.path.join(fmri_folder_path, fmri_subj)
        fmri_df = load_fmri_data_from_excel(fmri_subj_path)
        clean_fmri_df = clean_from_irrelevant_graph_metrics(fmri_df)

        # extract name of subj
        subj_name = fmri_subj.split('_')[0]

        combined_fmri_data['Net Nodal Strength'] = pd.concat([combined_fmri_data['Net Nodal Strength'], clean_fmri_df['Net NodStr'].rename(subj_name + "_Net NodStr")], axis=1)
        combined_fmri_data['Norm Nodal Strength'] = pd.concat([combined_fmri_data['Norm Nodal Strength'], clean_fmri_df['Norm NodStr'].rename(subj_name + "_Norm NodStr")], axis=1)
        combined_fmri_data['Positive Assortativity'] = pd.concat([combined_fmri_data['Positive Assortativity'], clean_fmri_df['Pos Assort'].rename(subj_name + "_Pos Assort")], axis=1)
        combined_fmri_data['Local Efficiency'] = pd.concat([combined_fmri_data['Local Efficiency'], clean_fmri_df['Loc Eff'].rename(subj_name + "_Loc Eff")], axis=1)
        combined_fmri_data['Nodal Path Length'] = pd.concat([combined_fmri_data['Nodal Path Length'], clean_fmri_df['Nod PathLen'].rename(subj_name + "_Nod PathLen")], axis=1)
        combined_fmri_data['Nodal Betweenness'] = pd.concat([combined_fmri_data['Nodal Betweenness'], clean_fmri_df['Nod Betweenness'].rename(subj_name + "_Nod Betweenness")], axis=1)

    return combined_fmri_data


results_master = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\2_CrossModalityGroupAnalysis\3_NodalMetricComparison"
if not os.path.exists(results_master):
    os.makedirs(results_master)

# folder paths
fmri_biv_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\6_GraphMetricsPerSubject\1_StandardCor"
fnirs_biv_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\4_GraphMetricsPerSubject\1_StandardCor"


fmri_par_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\6_GraphMetricsPerSubject\2_PartialCor"
fnirs_par_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\4_GraphMetricsPerSubject\2_PartialCor"

# load files
fmri_biv_files = [f for f in os.listdir(fmri_biv_path) if f.endswith('.xlsx')]
fnirs_biv_files = [f for f in os.listdir(fnirs_biv_path) if f.endswith('.xlsx')]
fmri_par_files = [f for f in os.listdir(fmri_par_path) if f.endswith('.xlsx')]
fnirs_par_files = [f for f in os.listdir(fnirs_par_path) if f.endswith('.xlsx')]

# Now we want to extract the node i from each file and perform ttest across subjects

# Bivariate First
# combine data into dictionaries
print("Loading Bivariate Data...")
fnirs_biv_hbo_dict, fnirs_biv_hbr_dict = combine_fnirs_data_to_dict(fnirs_biv_path)

fmri_biv_dict = combine_fmri_data_to_dict(fmri_biv_path)

print("Bivariate Analysis")
print("fMRI vs fNIRS HbO")
txt.append("Bivariate Analysis fMRI vs fNIRS HbO")

biv_hbo_ttest_dict = perform_ttest_dictionary(fnirs_biv_hbo_dict, fmri_biv_dict)
print("fMRI vs fNIRS HbR")
txt.append("Bivariate Analysis fMRI vs fNIRS HbR")

biv_hbr_ttest_dict = perform_ttest_dictionary(fnirs_biv_hbr_dict, fmri_biv_dict)
# save results to excel and each metric on different sheet
biv_hbo_out = os.path.join(results_master, "GraphMetric_nodal_ttest_fmri_vs_fnirs_bivariate_hbo.xlsx")
with pd.ExcelWriter(biv_hbo_out) as writer:
    for metric, result_df in biv_hbo_ttest_dict.items():
        result_df.to_excel(writer, sheet_name=metric)
biv_hbr_out = os.path.join(results_master, "GraphMetric_nodal_ttest_fmri_vs_fnirs_bivariate_hbr.xlsx")
with pd.ExcelWriter(biv_hbr_out) as writer:
    for metric, result_df in biv_hbr_ttest_dict.items():
        result_df.to_excel(writer, sheet_name=metric)
print("Bivariate t-test results saved.")

# Partial Correlation Next
print("Loading Partial Correlation Data...")
fnirs_par_hbo_dict, fnirs_par_hbr_dict = combine_fnirs_data_to_dict(fnirs_par_path)
fmri_par_dict = combine_fmri_data_to_dict(fmri_par_path)

print("Partial Correlation Analysis")
print("fMRI vs fNIRS HbO")
txt.append("Partial Correlation Analysis fMRI vs fNIRS HbO")
par_hbo_ttest_dict = perform_ttest_dictionary(fnirs_par_hbo_dict, fmri_par_dict)
print("fMRI vs fNIRS HbR")
txt.append("Partial Correlation Analysis fMRI vs fNIRS HbR")
par_hbr_ttest_dict = perform_ttest_dictionary(fnirs_par_hbr_dict, fmri_par_dict)
# save results to excel and each metric on different sheet
par_hbo_out = os.path.join(results_master, "GraphMetric_nodal_ttest_fmri_vs_fnirs_partialcor_hbo.xlsx")
with pd.ExcelWriter(par_hbo_out) as writer:
    for metric, result_df in par_hbo_ttest_dict.items():
        result_df.to_excel(writer, sheet_name=metric)   
par_hbr_out = os.path.join(results_master, "GraphMetric_nodal_ttest_fmri_vs_fnirs_partialcor_hbr.xlsx")
with pd.ExcelWriter(par_hbr_out) as writer:
    for metric, result_df in par_hbr_ttest_dict.items():
        result_df.to_excel(writer, sheet_name=metric)

print("*"*50)
for line in txt:
    print(line)

# save the text in the output folder
summary_txt_path = os.path.join(results_master, "ttest_summary.txt")
with open(summary_txt_path, 'w') as f:
    for line in txt:
        f.write(line + '\n')