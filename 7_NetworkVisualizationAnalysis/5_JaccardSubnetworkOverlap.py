"""
In this script we are going to check the overlap between fnirs and fmri modules with set theory (Jacard index)


for example:
for fmri module A, how many nodes are also in fnirs module A, B, C, etc
take the maximum overlap as the matching module
save the results in a text or excel
"""

import os
import pandas as pd
import numpy as np
from sklearn.metrics import jaccard_score
import time
def remove_irrelevant_columns(df):
    # we only keep the columns that contain both 'ROI' and 'UC=7' (unique comminities = 7)
    df = df.copy()
    for col in df.columns:
        if 'Pair Satori' not in col and 'UC=7' not in col:
            df = df.drop(col, axis=1, errors='ignore')
    return df

def community_dict_from_df(df):
    # The community labels are in the second column
    df = df.copy()
    community_col = df.columns[1]
    roi_col = df.columns[0]
    comm_dict = {}
    for _, row in df.iterrows():
        comm = row[community_col]
        roi = row[roi_col]
        # print(comm, roi)
        # quit()
        if comm not in comm_dict:
            comm_dict[comm] = set()
        comm_dict[comm].add(roi)

    return comm_dict

def jaccard_calculation(fmri_dict, fnirs_dict):

    results = []
    for fmri_comm, fmri_nodes in fmri_dict.items():

        # print(f"FMRI Community {fmri_comm}:")
        # print(f"    Nodes: {fmri_nodes}")
        # print(f"checking overlap of fMRI community {fmri_comm} with fNIRS HbO communities:")
        
        for fnirs_comm, fnirs_nodes in fnirs_dict.items():
            # print(f"    fNIRS HbO Community {fnirs_comm}:")
            # print(f"            Nodes: {fnirs_nodes}")

            intersection = fmri_nodes.intersection(fnirs_nodes)
            # print(f"            Overlapping Nodes: {intersection}")
            union = fmri_nodes.union(fnirs_nodes)
            jaccard_index = len(intersection) / len(union) if len(union) > 0 else 0 # manual calc of jaccard, because why not
            # print(f"   Jaccard Index: {jaccard_index:.3f}")

            results.append({
                "fMRI Community Label": fmri_comm,
                "fNIRS Community Label": fnirs_comm,
                "JaccardIndex": jaccard_index
            })
    # Convert the results list to a DataFrame
    results_df = pd.DataFrame(results)
    return results_df      

# sort it by fMRI comm and then by jaccard
def find_max_jaccard(fmri_fnirs_jac_df):
    df = fmri_fnirs_jac_df.copy()
    sorted_df = df.sort_values(
        by=['fMRI Community Label', 'JaccardIndex'], ascending=[True, False]
        )
    # Group by fMRI Community Label and find the max JaccardIndex
    max_indices = sorted_df.groupby('fMRI Community Label')['JaccardIndex'].idxmax()
    max_sorted_df = sorted_df.loc[max_indices]

    # Safety check: Ensure no duplicate max JaccardIndex for any fMRI Community Label
    for comm_label, group in sorted_df.groupby('fMRI Community Label'):
        max_jaccard = group['JaccardIndex'].max()

        # if (group['JaccardIndex'] == max_jaccard).sum() > 1:
        #     print(f"Warning: Multiple max Jaccard indices found for fMRI Community '{comm_label}'")
        #     print(group[group['JaccardIndex'] == max_jaccard])

    return max_sorted_df.reset_index(drop=True)

def permutation_test_max_jaccard(fmri_dict, fnirs_dict_hbo, fnirs_dict_hbr, max_jaccard_fmri_hbo, max_jaccard_fmri_hbr):
    """
    In this function we will perform a permutation test to assess the significance of the observed maximum Jaccard indices between fMRI and fNIRS communities.
    By keeping fmri community assignments fixed and randomly shuffling fNIRS community assignments, we can create a null distribution of maximum Jaccard indices.
    and test the hypothesis that the observed overlaps of fNIRS with fMRI are greater than would be expected by chance.

    Steps:
    1. get the list of fMRI community labels with their ROIs
    2. for each permutation
        a. shuffle the fnirs ROIs among the fnirs communities
        b. recalculate the jaccard indices between fmri and shuffled fnirs communities
        c. find the max jaccard for each fmri community
        d. store the max jaccard indices for this permutation
    3. after all permutations, compare the observed max jaccard indices to the null distribution to compute p-values
    """
    # make the fmri dictionary into a list: fmri_list = [(comm, nodes), ...]
    fmri_list = list(fmri_dict.items())

    # for fnirs, we need to do something a bit different...\
    # lets extract all the fNIRS ROIs into a single list
    all_fnirs_rois = []
    for nodes in fnirs_dict_hbo.values():
        all_fnirs_rois.extend(list(nodes))

    # so basically we have all the fnirs ROIs in a single list, and so we can assign them randomly to the fnirs communities
    fnirs_communities = list(fnirs_dict_hbo.keys())

    # so in each itteration, we will randomly assign the ROIs to the communities and at the end see if the max jaccard is higher than the observed one. this can be done once, it doesnt matter between hbo and hbr, because we got a random assignment
    print("Starting permutation test...")
    n_iterations = 100000
    null_per_community = {comm: [] for comm in fmri_dict.keys()}
    # print(null_per_community)

    start = time.time()
    for i in range(n_iterations):
        # shuffle the fnirs ROIs among the fnirs communities
        # to be true, we need to shuffle each time the original and preserve its structure, otherwise permutation is not correct
        shuffled_rois = all_fnirs_rois.copy()
        np.random.shuffle(shuffled_rois)
        # now we need to assign them to the fnirs communities in a random manner, but preserving the length of the original communities 
        shuffled_fnirs_dict = {}
        start_idx = 0
        for comm in fnirs_communities:
            comm_size = len(fnirs_dict_hbo[comm])
            shuffled_fnirs_dict[comm] = set(shuffled_rois[start_idx:start_idx + comm_size])
            start_idx += comm_size
        # now we can calculate the jaccard indices between fmri and shuffled fnirs communities
        shuffled_jaccard_df = jaccard_calculation(fmri_dict=fmri_dict, fnirs_dict=shuffled_fnirs_dict)
        # find the max jaccard for each fmri community
        shuffled_max_jaccard_df = find_max_jaccard(shuffled_jaccard_df)

        # now we can build the null distribution of null max jaccard indices per community
        for _, row in shuffled_max_jaccard_df.iterrows():
            comm_label = row['fMRI Community Label']
            jaccard_index = row['JaccardIndex']
            null_per_community[comm_label].append(jaccard_index)
    end = time.time()
    elapsed = (end - start)/60
    print(f"Permutation test completed in {elapsed:.2f} minutes.")

    # theoretically we have our null distribution for fnirs, and its irrelevant if its hbo or hbr.
    # lets calculate the p-values for hbo and hbr observed max jaccard indices
    print("Calculating p-values for fNIRS HbO:")
    p_values_hbo = []
    for _, row in max_jaccard_fmri_hbo.iterrows():
        comm_label = row['fMRI Community Label']
        observed_jaccard = row['JaccardIndex']

        null_distribution = null_per_community[comm_label]

        p_greater = (np.sum(np.array(null_distribution) >= observed_jaccard) + 1) / (n_iterations + 1) # one sided p-value: greater overlap than chance level

        p_values_hbo.append({'fMRI Community Label': comm_label,
                            'Observed Jaccard': observed_jaccard,
                            'p_greater': p_greater,
                            })
    p_values_hbo = pd.DataFrame(p_values_hbo)
    print(p_values_hbo)

    print("Calculating p-values for fNIRS HbR:")
    p_values_hbr = []
    for _, row in max_jaccard_fmri_hbr.iterrows():
        comm_label = row['fMRI Community Label']
        observed_jaccard = row['JaccardIndex']

        null_distribution = null_per_community[comm_label]

        p_greater = (np.sum(np.array(null_distribution) >= observed_jaccard) + 1) / (n_iterations + 1) # one sided p-value: greater overlap than chance level

        p_values_hbr.append({'fMRI Community Label': comm_label,
                            'Observed Jaccard': observed_jaccard,
                            'p_greater': p_greater,
                            })
    p_values_hbr = pd.DataFrame(p_values_hbr)
    print(p_values_hbr)

    return p_values_hbo, p_values_hbr

# load the fnirs and fmri module assignment files
out= r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\2_CrossModalityGroupAnalysis\4_ModularityOverlaps"
if not os.path.exists(out):
    os.makedirs(out)

fmri_files_path =r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\7_GraphMetricsGroupLevel_Updated\BrainNetFiles\1_RawDataFrames"
fnirs_files_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\5_GraphMetricsGroupLevel_Updated\BrainNetFiles\1_RawDataFrames"

fmri_excel_files = [f for f in os.listdir(fmri_files_path) if f.endswith('DF.xlsx')]
fnirs_excel_files = [f for f in os.listdir(fnirs_files_path) if f.endswith('DF.xlsx')]

for fmri_file in fmri_excel_files:
    for fnirs_file in fnirs_excel_files:
        #find the matching files
        if "Bivariate" in fmri_file and "Bivariate" in fnirs_file:
            # print(f"Comparing {fmri_file} and {fnirs_file}")
            fmri_df = pd.read_excel(os.path.join(fmri_files_path, fmri_file))
            fnirs_oxy_df = pd.read_excel(os.path.join(fnirs_files_path, fnirs_file), sheet_name='HbO_Bivariate_Corr')
            fnirs_deoxy_df = pd.read_excel(os.path.join(fnirs_files_path, fnirs_file), sheet_name='HbR_Bivariate_Corr')
        elif "Partial" in fmri_file and "Partial" in fnirs_file:
            # print(f"Comparing {fmri_file} and {fnirs_file}")
            fmri_df = pd.read_excel(os.path.join(fmri_files_path, fmri_file))
            fnirs_oxy_df = pd.read_excel(os.path.join(fnirs_files_path, fnirs_file), sheet_name='HbO_Partial_Corr')
            fnirs_deoxy_df = pd.read_excel(os.path.join(fnirs_files_path, fnirs_file), sheet_name='HbR_Partial_Corr')
        else:
            print(f"Skipping non-matching files: {fmri_file} and {fnirs_file}")
            continue
        print(f"Comparing {fmri_file} and {fnirs_file}")
        # now the data are loaded and we can check the overlap
        # 1st we remove irrelevant columns; for testing we are going to use only UC_7
        print(fmri_df.columns)

        fmri_df = remove_irrelevant_columns(fmri_df)
        fnirs_oxy_df = remove_irrelevant_columns(fnirs_oxy_df)
        fnirs_deoxy_df = remove_irrelevant_columns(fnirs_deoxy_df)

        # now we can check the overlap of the communities
        # we will construct a set for each community
        # example
        # fmri_dict = {
        #     'A': set of nodes in community A,
        #     'B': set of nodes in community B,
        #     ...
        # }
        fmri_dict = community_dict_from_df(fmri_df)

        fnirs_oxy_dict = community_dict_from_df(fnirs_oxy_df)
        fnirs_deoxy_dict = community_dict_from_df(fnirs_deoxy_df)
        
        fmri_fnirs_hbo_jaccard_df = jaccard_calculation(fmri_dict=fmri_dict, fnirs_dict=fnirs_oxy_dict)
        fmri_fnirs_hbr_jaccard_df = jaccard_calculation(fmri_dict=fmri_dict, fnirs_dict=fnirs_deoxy_dict)

        #* Sanity ckeck: is jaccard correctly calculated?
        # sanity=jaccard_calculation(fmri_dict=fmri_dict, fnirs_dict=fmri_dict)
        # print("Sanity check (fMRI vs fMRI):")
        # print(sanity)
        # quit()
        #NOTE its correctly calculated

        max_jaccard_fmri_hbo = find_max_jaccard(fmri_fnirs_hbo_jaccard_df)
        max_jaccard_fmri_hbr = find_max_jaccard(fmri_fnirs_hbr_jaccard_df)

        # compute the total overlap score
        total_overlap_hbo = max_jaccard_fmri_hbo['JaccardIndex'].sum()/len(max_jaccard_fmri_hbo['JaccardIndex'])
        total_overlap_hbr = max_jaccard_fmri_hbr['JaccardIndex'].sum()/len(max_jaccard_fmri_hbr['JaccardIndex'])

        # PERMUTATION TEST COULD BE ADDED HERE
        # we need to randomly shuffle the fnirs commu
        p_hbo, p_hbr = permutation_test_max_jaccard(
            fmri_dict=fmri_dict, 
            fnirs_dict_hbo=fnirs_oxy_dict,
            fnirs_dict_hbr=fnirs_deoxy_dict,
            max_jaccard_fmri_hbo=max_jaccard_fmri_hbo,
            max_jaccard_fmri_hbr=max_jaccard_fmri_hbr
        )

        # add a row at the end of the df with the total overlap score
        total_row_hbo = pd.DataFrame({
            "fMRI Community Label": ["Total Overlap Score"],
            "fNIRS Community Label": [""],
            "JaccardIndex": [total_overlap_hbo]
        })
        total_row_hbr = pd.DataFrame({
            "fMRI Community Label": ["Total Overlap Score"],
            "fNIRS Community Label": [""],
            "JaccardIndex": [total_overlap_hbr]
        })
        max_jaccard_fmri_hbo = pd.concat([max_jaccard_fmri_hbo, total_row_hbo], ignore_index=True)
        max_jaccard_fmri_hbr = pd.concat([max_jaccard_fmri_hbr, total_row_hbr], ignore_index=True)

        # append the p-values to the df
        # Merge the p-values (p_hbo) with the max_jaccard_fmri_hbo DataFrame
        max_jaccard_fmri_hbo = max_jaccard_fmri_hbo.merge(
            p_hbo[['fMRI Community Label', 'p_greater']],  # Select only the relevant columns
            on='fMRI Community Label',  # Merge on the community label
            how='left'  # Use a left join to preserve all rows in max_jaccard_fmri_hbo
        )

        max_jaccard_fmri_hbr = max_jaccard_fmri_hbr.merge(
            p_hbr[['fMRI Community Label', 'p_greater']],  # Select only the relevant columns
            on='fMRI Community Label',  # Merge on the community label
            how='left'  # Use a left join to preserve all rows in max_jaccard_fmri_hbr
        )

        #! here is the remapping, if not needed, comment it out
        # Remap the community labels from letters to numbers for easier reading
        label_dict = {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7}
        max_jaccard_fmri_hbo['fMRI Community Label'] = max_jaccard_fmri_hbo['fMRI Community Label'].map(label_dict).fillna(max_jaccard_fmri_hbo['fMRI Community Label'])
        max_jaccard_fmri_hbo['fNIRS Community Label'] = max_jaccard_fmri_hbo['fNIRS Community Label'].map(label_dict).fillna(max_jaccard_fmri_hbo['fNIRS Community Label'])
        max_jaccard_fmri_hbr['fMRI Community Label'] = max_jaccard_fmri_hbr['fMRI Community Label'].map(label_dict).fillna(max_jaccard_fmri_hbr['fMRI Community Label'])
        max_jaccard_fmri_hbr['fNIRS Community Label'] = max_jaccard_fmri_hbr['fNIRS Community Label'].map(label_dict).fillna(max_jaccard_fmri_hbr['fNIRS Community Label'])

        # save the results
        if "Bivariate" in fmri_file:
            con_type = "Bivariate"
        elif "Partial" in fmri_file:
            con_type = "Partial"
        output_hbo_file = os.path.join(out, f"fmri_fnirs_hbo_jaccard_{con_type}.xlsx")
        output_hbr_file = os.path.join(out, f"fmri_fnirs_hbr_jaccard_{con_type}.xlsx")
        max_jaccard_fmri_hbo.to_excel(output_hbo_file, index=False)
        max_jaccard_fmri_hbr.to_excel(output_hbr_file, index=False)

print("Results saved in:", out)




       
               


            




