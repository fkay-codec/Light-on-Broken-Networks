"""
In this script we will compare the averaged correlation matrices between fNIRS and fMRI
using the Mantel test.

How the Mantel test works:
-	You have two matrices (say, fNIRS FC and fMRI FC).
Step 1: Flatten the upper triangles into vectors.
Step 2: Compute the correlation between the two vectors (Pearson or Spearman). Call this the observed statistics
Step 3: Permutation step:
-	Shuffle the row/column labels of one matrix many times (e.g., 10,000 permutations).
-	Each time, recompute the correlation between the permuted matrix and the other matrix.
-	This builds a null distribution of correlations that would be expected by chance if the matrices were unrelated.
-	P-value: is the probability of obtaining a correlation at least as strong as your observed one, under the null hypothesis that the two matrices are unrelated
-	P<0.05: significantly more similar than chance	
Its more valid than anything else because it accounts for the fact that the matrix variables are not independent 

INPUT:
    - Ztransformed averaged correlation matrices from fNIRS and fMRI (output of 2_AverageCorM.py in both modalities)

OUTPUT:
    - Text file with the Mantel test results (correlation coefficient and p-value)
#! Verified that it is correct with 2a_spearans.py || it is basically the same but the permutations are done automaticlally and there is no need to do the upper yourself... 
#! Verified everything is correct... :D
"""
import os
import pandas as pd
import scipy.stats as stats
from skbio import DistanceMatrix
from skbio.stats.distance import mantel
from scipy.spatial.distance import squareform
import numpy as np

def isMatch(df1, df2):
    """Check if two dataframes have the same shape and labels"""
    if df1.shape != df2.shape:
        return False
    if not df1.columns.equals(df2.columns):
        return False
    if not df1.index.equals(df2.index):
        return False
    return True

def rename_fnirs_df(fnirs_df):
    """replace '_' with '' in both columns and index of the dataframe"""
    fnirs_df.columns = [col.replace ("_", "") for col in fnirs_df.columns]
    fnirs_df.index = [idx.replace ("_", "") for idx in fnirs_df.index]

    return fnirs_df

def rename_fmri_df(fmri_df):
    """"
    Replace the ROI1_ from the columns and index of the dataframe
    """
    fmri_df.columns = [col.split("_")[-1] for col in fmri_df.columns]
    fmri_df.index = [idx.split("_")[-1] for idx in fmri_df.index]

    return fmri_df

def inverse_z_transform(inverse_z_df):
    """Inverse Fisher Z-transformation from the dataframe"""
    # make the diagonal 1
    inversed_z = np.tanh(inverse_z_df)
    np.fill_diagonal(inversed_z.values, 1)
    return inversed_z

def similarity_to_distance(df):
    """
    Convert a similarity matrix (correlation) to a distance matrix
    We want to create distance matrices because the mantel function from skbio
    requires distance matrices as input.

    In our analysis we have computed the correlation matrices which are essentially similarity matrices.

    To convert them to distance we need to do: distance = 1 - similarity
    Note: in this case we should ask ourselves, how should we treat negative correlations?
        Either we treat them as equally distant with the positives
            eg. -0.8 and 0.8 are equally distant from 0
        Or we treat them as more distant than the positives
            eg. -0.8 is more distant from 0 than 0.8
        
        If we treated them as equally distance, in the fMRI dataset a -0.9 correlation will be treated as 1.9 distance and in the fNIRS the 0.9 correlation as 1.9 as well.
        This means that when we perform the Mantel's test, we have to consider what to do...

        In our case, evaluating functional connectivity, a negative correlation in fMRI and a positive in fNIRS indicates anti-synchronous activity... 
        If we keep the sign as is, we will asnwer the hypothesis:
            Are the connectivity patterns in fNIRS and fMRI similar, including the directionality of the connections?
        
    
    What we do:
        distance = 1 - similarity
            it is more valid for our case  
    """

    return 1-df


def upper(df):
    """Return upper triangle values (k=1) from DataFrame or ndarray."""
    try:
        assert(type(df) == np.ndarray)
    except:
        if type(df) == pd.DataFrame:
            df = df.values
        else:
            raise TypeError('Must be np.ndarray or pd.DataFrame')
    mask = np.triu_indices(df.shape[0], k=1)
    return df[mask]

# Define input and output paths
output_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\2_CrossModalityGroupAnalysis\1_MantelTestResults"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

input_fmri_master_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\4_AverageCorM"
input_fnirs_master_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\2_AverageCorM"


txt_results = []
txt_results.append("="*100 + "\n")
txt_results.append("Assessing similarity between fNIRS and fMRI correlation matrices using Mantel test\n")
txt_results.append("="*100 + "\n")

# Iterate through files in the fMRI folder
for fmri_file in os.listdir(input_fmri_master_folder):
    if not fmri_file.endswith(".xlsx"):  # Skip non-Excel files
        continue

    # Build the full path for the fMRI file
    fmri_file_path = os.path.join(input_fmri_master_folder, fmri_file)

    # Build the corresponding fNIRS file path
    # they have the same name so it is easy to match them
    fnirs_file_path = os.path.join(input_fnirs_master_folder, fmri_file)

    # Check if the corresponding fNIRS file exists
    if not os.path.exists(fnirs_file_path):
        print(f"    {fmri_file}: Corresponding fNIRS file not found. Exiting script...")
        quit()

    print(f"Processing Mantel test for: \n{fmri_file_path} \n{fnirs_file_path}\n...")

    txt_results.append(f"Mantel test for: \n\tfMRI: {os.path.basename(fmri_file_path).split('_')[-2]} Vs fNIRS: {os.path.basename(fnirs_file_path).split('_')[-2]}\n")


    # Load the fMRI and fNIRS data
    fmri_df = pd.read_excel(fmri_file_path, index_col=0)
    fnirs_oxy_df = pd.read_excel(fnirs_file_path, index_col=0, sheet_name="HbO_Average")
    fnirs_deoxy_df = pd.read_excel(fnirs_file_path, index_col=0, sheet_name="HbR_Average")

    # Rename columns and indices to ensure compatibility
    fmri_df = rename_fmri_df(fmri_df)
    fnirs_oxy_df = rename_fnirs_df(fnirs_oxy_df)
    fnirs_deoxy_df = rename_fnirs_df(fnirs_deoxy_df)
    # Sanity check: Ensure that the matrices are equivelent in terms of labels/ordering/size
    if not isMatch(fmri_df, fnirs_oxy_df) or not isMatch(fmri_df, fnirs_deoxy_df):
        print("The matrices do not match, and they SHOULD exiting script...")
        print(fmri_df.head())
        print(fnirs_oxy_df.head())
        print(fnirs_deoxy_df.head())
        quit()
    print("everything good, proceeding...")



    # --- Standard Pearson + permutation testing (save separately) ---
    output_folder_standard = os.path.join(output_folder, "StandardCorrPermutation")
    if not os.path.exists(output_folder_standard):
        os.makedirs(output_folder_standard)
    # Use Fisher z-values for Pearson + permutation — DO NOT inverse here
    fmri_upper = upper(fmri_df)
    fnirs_oxy_upper = upper(fnirs_oxy_df)
    fnirs_deoxy_upper = upper(fnirs_deoxy_df)

    # analytic Pearson
    cor_oxy, p_analytic_oxy = stats.pearsonr(fmri_upper, fnirs_oxy_upper)
    cor_deoxy, p_analytic_deoxy = stats.pearsonr(fmri_upper, fnirs_deoxy_upper)

    # permutation test: shuffle matrix labels (preserves matrix dependency)
    np.random.seed(0)
    n_iter = 100000
    m_ids = list(fmri_df.columns)
    m2_oxy = fnirs_oxy_upper.copy()
    m2_deoxy = fnirs_deoxy_upper.copy()
    rhos_oxy = []
    rhos_deoxy = []
    for i in range(n_iter):
        np.random.shuffle(m_ids)
        perm_fmri_upper = upper(fmri_df.loc[m_ids, m_ids])
        r_oxy = stats.pearsonr(perm_fmri_upper, m2_oxy)[0]
        r_deoxy = stats.pearsonr(perm_fmri_upper, m2_deoxy)[0]
        rhos_oxy.append(r_oxy)
        rhos_deoxy.append(r_deoxy)

    perm_p_oxy = ((np.sum(np.abs(cor_oxy) <= np.abs(rhos_oxy)))+1)/(n_iter+1)
    perm_p_deoxy = ((np.sum(np.abs(cor_deoxy) <= np.abs(rhos_deoxy)))+1)/(n_iter+1)

    print("Standard Pearson + permutation results:")
    print(f"  HbO: r = {cor_oxy:.4e}, analytic p = {p_analytic_oxy:.2e}, perm p = {perm_p_oxy:.2e}")
    print(f"  HbR: r = {cor_deoxy:.4e}, analytic p = {p_analytic_deoxy:.2e}, perm p = {perm_p_deoxy:.2e}")

    # save to separate file
    standard_txt = []
    standard_txt.append(f"Standard Pearson + permutation for: {os.path.basename(fmri_file)}\n")
    standard_txt.append(f"HbO: r = {round(cor_oxy, 4)}, analytic p = {p_analytic_oxy:.2e}, perm p = {perm_p_oxy:.2e}\n")
    standard_txt.append(f"HbR: r = {round(cor_deoxy, 4)}, analytic p = {p_analytic_deoxy:.2e}, perm p = {perm_p_deoxy:.2e}\n")
    standard_txt.append(f"Permutations: {n_iter}\n")

    with open(os.path.join(output_folder_standard, f"StandardCorr_{os.path.basename(fmri_file)}.txt"), "w") as f:
        f.writelines(standard_txt)
    # ---------------------------------------------------------------
    # Prepare the data for the Mantel test
    # inverse Fisher Z-transformation to get back to correlation values
    fmri_df = inverse_z_transform(fmri_df)
    fnirs_oxy_df = inverse_z_transform(fnirs_oxy_df)
    fnirs_deoxy_df = inverse_z_transform(fnirs_deoxy_df)
    # convert the similarity matrices to distance matrices
    fmri_df = similarity_to_distance(fmri_df)
    fnirs_oxy_df = similarity_to_distance(fnirs_oxy_df)
    fnirs_deoxy_df = similarity_to_distance(fnirs_deoxy_df)

    # fmri_dm = DistanceMatrix(fmri_df.values, ids=fmri_df.index)
    # fnirs_oxy_dm = DistanceMatrix(fnirs_oxy_df.values, ids=fnirs_oxy_df.index)
    # fnirs_deoxy_dm = DistanceMatrix(fnirs_deoxy_df.values, ids=fnirs_deoxy_df.index)
    # print(fmri_dm)

    print("Running Mantel test for HbO...")
    corr1, p_value1, _ = mantel(fmri_df.values, fnirs_oxy_df.values, method="pearson", permutations=100000, seed=42)
    print(f"HbO Mantel test results: correlation={round(corr1, 4)}, p-value={p_value1:.2e}")
    txt_results.append(f"Similarity between fMRI and fNIRS HbO correlation matrices:\nr = {round(corr1, 4)}, p = {p_value1:.2e}\n")
    

    print("Running Mantel test for HbR...")
    corr2, p_value2, _ = mantel(fmri_df.values, fnirs_deoxy_df.values, method="pearson", permutations=100000, seed=42)
    print(f"HbR Mantel test results: correlation={round(corr2, 4)}, p-value={p_value2:.2e}")
    txt_results.append(f"Similarity between fMRI and fNIRS HbR correlation matrices:\nr = {round(corr2, 4)}, p = {p_value2:.2e}\n")
    txt_results.append("-"*50 + "\n\n")



txt_results.append("End of the results\n")
txt_results.append("Mantel's parameters:\n\tMethod: Pearson\n\tPermutations: 100000\n\tRandom seed: 42\n")
with open(os.path.join(output_folder, "MantelTest_Results.txt"), "w") as f:
    f.writelines(txt_results)

print("All done!")