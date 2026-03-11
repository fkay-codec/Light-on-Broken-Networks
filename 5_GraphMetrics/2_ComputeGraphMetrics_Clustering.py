"""
================================================================================
Script: Group-Level Graph Metrics and Modularity Computation
================================================================================

PURPOSE
--------
This script computes graph-theoretical metrics and community structure (modularity) 
for group-level correlation matrices derived from fMRI and fNIRS data. It processes 
Fisher z-transformed correlation matrices, computes nodal-level graph metrics, and 
identifies stable community assignments across multiple resolution parameters (gamma). 
-----------------------------------------------------------------------------------
Input Data:
    1. Group-Averaged Correlation Matrices:
       - fMRI: Located in `fmri_folder`
       - fNIRS: Located in `fnirs_folder` (separate sheets for HbO and HbR)
       - All matrices are Fisher Z-transformed and globally normalized to [0, 1] range
       - Normalization: divided by global maximum across ALL subjects and modalities
       - This enables direct cross-modal comparison of graph metrics

    2. Format:
       - Excel files (.xlsx) with group-averaged correlation matrices
       - Row/column indices represent brain regions or channels (N=82)
       - Diagonal values set to 0 before all computations 

================================================================================
GRAPH METRICS COMPUTATION
--------------------------
The script computes nodal-level graph metrics from group-averaged correlation matrices:

1. **Nodal Strength Variants**
   - Net, positive, negative, absolute, and normalized nodal strength
   - Number of positive and negative edges per node

2. **Local Assortativity**
   - Positive and negative assortativity coefficients
   - Computed using brainconn library (handles signed weights)

3. **Shortest Path Lengths**
   - Distance metric: 1 - |z| where z is Fisher Z-transformed correlation
   - Treats positive and negative correlations equally (both strong correlations = short distances)

4. **Betweenness Centrality**
   - Uses same distance metric as path length: 1 - |z|
   - Measures node importance based on shortest paths

5. **Local Efficiency**
   - Uses absolute values: |z| (treats negative correlations as positive)
   - BCT algorithm requires positive weights in [0, 1] range

Distance Metric Rationale:
   - Using 1 - |z| treats strong negative correlations (e.g., r=-0.9) as equally 
     "close" to strong positive correlations (e.g., r=0.9)
   - Both would have distance = 0.1, reflecting strong coupling regardless of sign
   - Alternative approaches (e.g., treating negatives as disconnected) were considered 
     but not implemented in this final version

================================================================================
MODULARITY COMPUTATION (CONSENSUS CLUSTERING APPROACH)
-------------------------------------------------------
Target: Find stable community structures with 5-7 communities (based on known 
cortical functional systems).

Gamma Parameter:
  - Controls community size (resolution parameter)
  - Gamma < 1 -> Larger communities
  - Gamma > 1 -> Smaller communities
  - Searched range: 0.3 to ~2.0 (terminates when consistently >10 communities)

CHALLENGE: Louvain Algorithm Stochasticity
-------------------------------------------
The Louvain algorithm is non-deterministic:
  - Each run can assign DIFFERENT numeric labels to the SAME community structure
  - Example:
      Run 1: Node_5 -> Community 2  
      Run 2: Node_5 -> Community 4  (same structural community, different label)
  - Single-run results are unreliable for determining stable community structure

SOLUTION: Consensus Clustering
-------------------------------
1. For each promising gamma value, run Louvain 1000 times (random seeds)
2. Build co-assignment matrix: tracks how often each node pair shares a community
3. Apply consensus algorithm (netneurotools.modularity.find_consensus):
   - Constructs agreement matrix from node pair co-occurrences
   - Clusters this agreement matrix to extract most stable community structure
4. Result: Single consensus community assignment representing 1000 stochastic runs

OUTPUT SELECTION STRATEGY:
--------------------------
For each target community count (5, 6, 7):
  - Store the gamma with HIGHEST median Q (modularity score)
  - Q measures within-community edge density above chance
  - Higher Q = stronger, more reliable community structure
  - Reject singleton communities (isolated single nodes)

THREE-STAGE COMPUTATIONAL EFFICIENCY PIPELINE:
-----------------------------------------------
Stage 1 (Quick Screen):
  - Run consensus with 50 iterations per gamma
  - Purpose: Rapidly identify gamma values yielding ~5-7 communities
  - Accepts range: 3-9 communities (loose filter)


Stage 2 (Verification):
  - Run consensus with 100 iterations
  - Purpose: Confirm gamma produces target range (5-7 communities)
  - Filters to exact target range


Stage 3 (Full Consensus):
  - Run 1000 Louvain iterations + consensus clustering
  - Purpose: Generate final stable community assignments
  - Tracks median Q for quality comparison
  - Only executed for gamma values passing Stage 2


Termination: When gamma consistently produces >10 communities (verified with 
100 iterations), search stops (further increases only fragment communities more).

"""


import warnings
import brainconn.core as bc

# Suppress RuntimeWarnings globally
warnings.filterwarnings("ignore", category=RuntimeWarning)
import string
import os
import pandas as pd
import numpy as np
import bct.algorithms as bct
from scipy.stats import mode
import time
from netneurotools import modularity
#! the difference with the other script lies in the search_best_gamma_seed_combination function in which the stable community the Q refreshes... see below at the #!
def search_best_gamma_seed_combination(data_array, gamma_step=0.05):
    """
    Find optimal gamma values yielding stable 5, 6, and 7-community structures
    using consensus clustering across 1000 Louvain iterations.
    
    Parameters:
    -----------
    data_array : np.ndarray
        Group-averaged correlation matrix (Fisher Z-transformed, globally normalized)
    gamma_step : float
        Gamma increment for search (default: 0.05) For final analysis we used 0.0125
    
    Returns:
    --------
    mod_dataframe : pd.DataFrame
        Three columns with consensus community assignments:
        - Column format: "Q={value}|G={gamma}|UC={count}"
        - One column each for 5, 6, and 7 communities
        - Values are relabeled community letters (A, B, C, ...)
    
    Algorithm (Three-Stage Pipeline):
    ---------------------------------
    
    Initialization:
        - Start at gamma = 0.3
        - Increment by gamma_step each iteration
    
    Stage 1 - Quick Screen (50 iterations):
        - Apply consensus clustering with 50 Louvain runs
        - If consensus yields 3-9 communities → proceed to Stage 2
        - If consensus yields >10 communities → verify termination
    
    Stage 2 - Verification (100 iterations):
        - Apply consensus clustering with 100 Louvain runs  
        - If consensus yields exactly 5, 6, or 7 communities -> proceed to Stage 3
    
    Stage 3 - Full Consensus (1000 iterations):
        - Run Louvain 1000 times, collect community assignments (ci) and Q scores
        - Apply consensus clustering to all 1000 assignments
        - Count unique communities in consensus result
        - Reject if any singleton communities (isolated nodes)
        - For each target count (5, 6, 7):
            * Store first valid result
            * Update if later gamma yields higher median Q
    
    Termination:
        - Stops when gamma consistently produces >10 communities (verified at 100 iters)
        - Rationale: Gamma >1 fragments communities; further increases won't help
    
    Output:
        - Returns best (highest median Q) consensus assignments for 5, 6, 7 communities
        - Raises ValueError if any target count not found
    """
    def verify_number_of_unique_communities(array, gamma, itterations):
        """
        Compute stable community count for given gamma via consensus clustering.
        
        Runs Louvain 'itterations' times, applies consensus algorithm to stabilize
        assignments, and returns the number of unique communities in consensus.
        """
        ci_list = []

        for i in range(itterations):
            ci, _ = bct.modularity_louvain_und_sign(array, gamma=gamma, seed=None)
            ci_list.append(ci)
        consensus = modularity.find_consensus(np.column_stack(ci_list))
        unique_consensus = np.unique(consensus).size

        return unique_consensus

    
    time_start = time.time()
    gamma = 0.3
    actual_unique_com = None
    mod_dataframe = pd.DataFrame()


    # we initialize the lists to store our arrays of CIs and Qs for each target community
    ci5_stable=None
    ci6_stable=None
    ci7_stable=None
    q5_stable=None
    q6_stable=None
    q7_stable=None




    while True:
        # Stage 1: Quick screen with 50-iteration consensus
        # Avoids expensive 1000-iteration consensus for unsuitable gamma values
        ci_quick = verify_number_of_unique_communities(data_array, gamma, itterations=50)

        # Termination check: if >10 communities in quick screen, verify with more iterations
        if ci_quick > 10:
            print("    Quick screen shows >10 communities, verifying...")
            ci_terminate = verify_number_of_unique_communities(data_array, gamma, itterations=100)
            if ci_terminate >= 11:
                print(f"    Confirmed: {ci_terminate} communities at gamma={gamma}. Terminating search.")
                break


        # Stage 2: If quick screen in range [3-9], verify with 100-iteration consensus
        if 3 <= ci_quick <= 9:
            # More stable estimate with 100 iterations
            ci_verify = verify_number_of_unique_communities(data_array, gamma, itterations=100)

            # Stage 3: If verified range [5-7], run full 1000-iteration consensus
            if 5 <= ci_verify <= 7:
                # Run 1000 Louvain iterations with random seeds
                ci_list = []
                q_list = []
                for i in range(1000):
                    actual_ci, actual_q = bct.modularity_louvain_und_sign(data_array, gamma=gamma, seed=None)
                    ci_list.append(actual_ci)
                    q_list.append(actual_q)
                
                # Apply consensus clustering to stabilize community assignments
                # Process:
                #   1. Build co-assignment matrix: tracks node pair co-occurrences across 1000 runs
                #   2. Cluster agreement matrix to extract most stable community structure
                #   3. Result: single consensus assignment representing all 1000 stochastic runs
                consensus = modularity.find_consensus(np.column_stack(ci_list))
                unique_cons = np.unique(consensus).size
                _, counts = np.unique(consensus, return_counts=True)

                # Reject solutions with singleton communities (isolated single nodes)
                if 1 in counts:
                    print(f"    Warning: Gamma {gamma} has singleton communities. Rejecting.")
                else:
                    # Storage logic: keep first valid result, update only if higher median Q found
                    if unique_cons == 5:
                        if ci5_stable is None:
                            print(f"    Storing CI for 5 communities at gamma {gamma}")
                            ci5_stable = consensus
                            q5_stable = np.median(q_list)
                            g5 = gamma
                        else:
                            if q5_stable < np.median(q_list):
                                print(f"    Updating CI for 5 communities at gamma {gamma} with higher Q; from {q5_stable} to {np.median(q_list)}")
                                ci5_stable = consensus
                                q5_stable = np.median(q_list)
                                g5 = gamma
                    if unique_cons == 6:
                        if ci6_stable is None:
                            print(f"    Storing CI for 6 communities at gamma {gamma}")
                            ci6_stable = consensus
                            q6_stable = np.median(q_list)
                            g6 = gamma
                        else:
                            if q6_stable < np.median(q_list):
                                print(f"    Updating CI for 6 communities at gamma {gamma} with higher Q; from {q6_stable} to {np.median(q_list)}")
                                ci6_stable = consensus
                                q6_stable = np.median(q_list)
                                g6 = gamma
                    
                    if unique_cons == 7:
                        if ci7_stable is None:
                            print(f"    Storing CI for 7 communities at gamma {gamma}")
                            ci7_stable = consensus
                            q7_stable = np.median(q_list)
                            g7 = gamma
                        else:
                            if q7_stable < np.median(q_list):
                                print(f"    Updating CI for 7 communities at gamma {gamma} with higher Q; from {q7_stable} to {np.median(q_list)}")
                                ci7_stable = consensus
                                q7_stable = np.median(q_list)
                                g7 = gamma


        print(f"Gamma: {gamma} processed.")
        gamma += gamma_step
        gamma= round(gamma, len(str(gamma_step).split(".")[-1])) # to avoid floating point issues; the second part basically converts the gamma step to string and counts the decimal points it has
    
    end_time = time.time()
    time_elapsed = end_time - time_start
    print(f"Search completed in {time_elapsed:.2f} seconds.")


    if all(x is not None for x in [q5_stable, g5, ci5_stable, q6_stable, g6, ci6_stable, q7_stable, g7, ci7_stable]):
        mod_dataframe[f"Q={q5_stable:.4f}|G={g5}|UC=5"] = ci5_stable
        mod_dataframe[f"Q={q6_stable:.4f}|G={g6}|UC=6"] = ci6_stable
        mod_dataframe[f"Q={q7_stable:.4f}|G={g7}|UC=7"] = ci7_stable
    else:
        raise ValueError("No valid gamma/seed combinations found for 5, 6, or 7 unique communities.\nEither broaden the verification criteria or lower the Gamma step")

    
    mod_dataframe_relab = relabel(mod_dataframe)




    return mod_dataframe_relab

def nodal_path_length(D):
    nodal_path_length = np.zeros(D.shape[0])

    # get the average distance while ignoring inf values (disconnected nodes), and self
    for i in range(D.shape[0]):
        mask = (np.arange(D.shape[0]) != i) & (D[i, :] != np.inf)
        nodal_path_length[i] = np.mean(D[i, mask])
    return nodal_path_length

def nodal_strength(df):
    df=df.copy()
    array = df.to_numpy()
    np.fill_diagonal(array, 0)
    df = pd.DataFrame(array, index=df.index, columns=df.columns)
    node_str = df.sum(axis=1) # all columns for each row
    return node_str

def inverse_z_transform(inverse_z_df):
    """Inverse Fisher Z-transformation from the numpy array"""
    # if the input is pandas dataframe convert to numpy array
    if isinstance(inverse_z_df, pd.DataFrame):
        inverse_z_df = inverse_z_df.to_numpy()
    # make the diagonal 1
    inversed_z = np.tanh(inverse_z_df)
    np.fill_diagonal(inversed_z, 1)
    return inversed_z

def fisher_z_transform(corr_matrix_df):
    """Apply Fisher z-transformation to a correlation matrix in Dataframe format"""

    clip = 0.9999999999  # To avoid infinities
    corr_matrix_df = corr_matrix_df.clip(-clip, clip)  # Avoid infinities
    z_matrix = np.arctanh(corr_matrix_df)
    return z_matrix

def modularity_array(corr_array, gamma, itterations=150):
    """
    Using the Louvain method (non deterministic) to compute modularity
    Modularity depends on gamma parameter || the default is 1 || large groups < 1 < small groups
    """
    
    ci_list = []
    q_list = []
    n_itter = itterations  # number of iterations to run the algorithm to get a stable result
    # corr_array = inverse_z_transform(corr_array)
    
    for i in range(n_itter):
        """
        example
        ci = [roi1: 1, roi2: 1, roi3: 2, roi4: 2, roi5: 3 ...] # community assignment
        """
        ci, q = bct.modularity_louvain_und_sign(
            corr_array, gamma=gamma, seed=42
        )
        ci_list.append(ci)
        q_list.append(q)
    # get the median modularity and the median community assignment
    q_median = np.median(q_list)
    # to find the most probable community assignment we are going to compute the mode
    # stack the list to a 2D pandas array [row-wise]
    ci_list = np.vstack(ci_list)
    df = pd.DataFrame(ci_list.T)
    outsave = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study"
    filepath = os.path.join(outsave, 'louvain_community_assignments_debug.xlsx')


    """
    Now that we have the list of community assignments for each iteration we have to group them based on their assignment
    the louvain algorithm each time can assign different community numbers to the same community.
    Eg. 
    - in iter 1 channel 1 was in commmunity 1
    - in iter 2 channel 1 was in community 3
    We have then to relabel them... the first community that occurs will be labeled A the second B etc. 
    This way we can then compute the mode across the rows to get the most probable community assignment
    """
    df2 = relabel(df)
    with pd.ExcelWriter(filepath) as writer:
        df.to_excel(writer, sheet_name='raw_out', index=False)
        df2.to_excel(writer, sheet_name='Relabeled_Communities', index=False)

    quit()

    df_mode = df.mode(axis=1)[0] # if they are ties we extract the first one || probably with 100 itterations there wont be ties || at 70, 100 and at 200 itters there are no ties (same as 1000) so n_itter = 150 better for safety

    return df_mode, q_median

def relabel(df):
    """
    Convert numeric community labels to letters for interpretability.

    Parameters:
    -----------
    df : pd.DataFrame
        Each column represents consensus community assignment from one gamma value
        Rows = nodes (channels), values = numeric community labels (1, 2, 3, ...)

    Returns:
    --------
    relabeled_df : pd.DataFrame
        Same structure as input, but numeric labels replaced with letters (A, B, C, ...)
    
    Purpose:
    --------
    Consensus clustering produces arbitrary numeric labels. Converting to letters:
    - Improves readability in output files
    - Makes community comparisons across gamma values easier
    - Standardizes labeling (first community = A, second = B, etc.)
    
    Method:
    -------
    For each column (gamma value):
      1. Identify unique numeric labels (e.g., [1, 3, 5])
      2. Map to letters in order of first appearance (1→A, 3→B, 5→C)
      3. Replace all instances in column
    """
    letters = list(string.ascii_uppercase + string.ascii_lowercase)  # Support up to 52 communities

    # Convert each column independently (each represents a different gamma)
    # Within each column: first unique label -> A, second -> B, etc.
    df_relabel = pd.DataFrame(index=df.index, columns=df.columns)
    for col in df.columns:
        unique_labels = df[col].unique()
        mapping = {old: letters[i] for i, old in enumerate(unique_labels)}
        df_relabel[col] = df[col].map(mapping)

    return df_relabel

def number_of_pos_neg_edges(corr_array):
    """
    Function to compute the number of positive and negative edges in a correlation matrix
    Input: correlation matrix in numpy array format
    Output: number of positive edges, number of negative edges
    """
    # Remove the diagonal and fill with 0
    np.fill_diagonal(corr_array, 0)

    # count the number of positives and negatives edge for each 'row'
    npos = np.sum(corr_array > 0, axis=1)
    nneg = np.sum(corr_array < 0, axis=1)

    return npos, nneg


def compute_graph_metrics(corr_df):
    """
    Compute node-level graph metrics from group-averaged correlation matrix.
    
    Input: 
        corr_df: DataFrame with Fisher Z-transformed, globally normalized correlation matrix [0, 1]
                 This is the group-averaged matrix (not individual subject)
    
    Output: 
        DataFrame with nodal metrics (one row per node/channel):
        - Net/Positive/Negative/Absolute/Normalized Nodal Strength
        - Number of positive/negative edges per node
        - Positive/Negative Local Assortativity
        - Local Efficiency using |z|
        - Nodal Path Length using distance 1-|z|
        - Betweenness Centrality using distance 1-|z|
    """

    corr_array = corr_df.to_numpy().copy()
    # Remove the diagonal and fill with 0
    np.fill_diagonal(corr_array, 0)


    # Compute nodal strength variants
    # Data is Fisher Z-transformed and globally normalized for cross-modal comparison
    # Global normalization preserves group mean differences (unlike within-matrix z-scoring)

    # Net nodal strength (sum of all connections)
    node_str = bct.strengths_und(corr_array)
    # Separate positive and negative contributions
    node_str_pos, node_str_neg, _, _ = bct.strengths_und_sign(corr_array)
    
    # Count number of positive and negative edges per node
    npos, nneg = number_of_pos_neg_edges(corr_array)
    # Sanity check: each node should have 81 edges total (82 nodes - 1 self)
    if np.any(npos+nneg != 81):
        print("Warning: number of positive and negative edges does not match the total number of edges")
        print("Exiting...")
        quit()
    
    # Normalized nodal strength based on Rubinov & Sporns (2011)
    # Accounts for balance between positive and negative connections
    # Dividing by 81 = (N-1) where N=82 nodes 
    norm_node_str = node_str_pos/81 - (node_str_neg/81)*(node_str_pos/(node_str_pos +node_str_neg))

    # Compute the absolute nodal strength
    node_str_abs = node_str_pos + np.abs(node_str_neg)

    # Compute local assortativity using brainconn (requires diagonal = 0)
    # Handles signed weights (positive and negative correlations)
    pos_assort, neg_assort = bc.local_assortativity_wu_sign(corr_array)


    # Compute shortest paths (requires converting correlation to distance matrix)
    # Matrices are in [0, 1] range after Fisher Z-transform and global normalization
    
    z = corr_df.to_numpy().copy()

    # Create distance matrix: L = 1 - |z|
    # Treats strong negative correlations as equally "close" to strong positives
    # (e.g., both r=0.9 and r=-0.9 → distance = 0.1)
    tmp = z.copy()
    tmp = np.abs(tmp)
    L1 = 1 - tmp

    # Compute all-pairs shortest path distances using BCT
    D1, _ = bct.distance_wei(L1)

    # Compute average shortest path length per node (excluding self and infinite distances)
    nodal_path_length1 = nodal_path_length(D1)

    # Betweenness centrality per node using distance metric 1-|z|
    # Measures node importance in shortest paths between other nodes
    nodal_betweeness1 = bct.betweenness_wei(L1)


    # Local efficiency using BCT algorithm
    # BCT requires positive weights in [0, 1], so we use absolute values
    # This treats negative correlations as positive with same magnitude
    tmp = z.copy()
    tmp = np.abs(tmp)
    local_eff = bct.efficiency_wei(tmp, local=True)

    # convert to pandas series for easier handling and saving
    # node str
    node_str = pd.Series(node_str, index=corr_df.index)
    node_str_pos = pd.Series(node_str_pos, index=corr_df.index)
    node_str_npos = pd.Series(npos, index=corr_df.index)
    node_str_neg = pd.Series(node_str_neg, index=corr_df.index)
    node_str_nneg = pd.Series(nneg, index=corr_df.index)
    node_str_abs = pd.Series(node_str_abs, index=corr_df.index)
    norm_node_str=pd.Series(norm_node_str, index=corr_df.index)
    # assortativity
    pos_assort = pd.Series(pos_assort, index=corr_df.index)
    neg_assort = pd.Series(neg_assort, index=corr_df.index)

    # path length
    nodal_path_length1 = pd.Series(nodal_path_length1, index=corr_df.index)


    # betweeness centrality
    nodal_betweeness1 = pd.Series(nodal_betweeness1, index=corr_df.index)
    
    # local efficiency
    local_eff= pd.Series(local_eff, index=corr_df.index)


    dataframe = pd.DataFrame({
        "Net Nodal Strength": node_str,
        "Positive Nodal Strength": node_str_pos,
        "Number of Positives": node_str_npos,
        "Negative Nodal Strength": node_str_neg,
        "Number of Negatives": node_str_nneg,
        "Absolute Nodal Strength": node_str_abs,
        "Normalized Nodal Strength": norm_node_str,
        "Positive Assortativity": pos_assort,
        "Negative Assortativity": neg_assort,
        "Local Efficiency [abs(z)]": local_eff,
        "Nodal Path Length (1-|z|)": nodal_path_length1,
        "Betweenness Centrality (1-|z|)": nodal_betweeness1,

    
    }, index=corr_df.index)    
    # Write the DataFrame to the Excel file
    return dataframe



# GROUP-LEVEL GRAPH METRICS AND MODULARITY
# Input: Group-averaged correlation matrices that are:
#   1. Fisher Z-transformed during averaging process
#   2. Globally normalized to [0, 1] by dividing by max across all subjects/modalities
#   3. This enables direct cross-modal comparison of graph metrics

fnirs_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\3_AverCorM_NormForGraph"
fnirs_output_folder = os.path.dirname(fnirs_folder) + r"\5_GraphMetricsGroupLevel_Updated"
if not os.path.exists(fnirs_output_folder):
    os.makedirs(fnirs_output_folder)
fmri_folder = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\5_AverCorM_NormForGraph"
fmri_output_folder = os.path.dirname(fmri_folder) + r"\7_GraphMetricsGroupLevel_Updated"
if not os.path.exists(fmri_output_folder):
    os.makedirs(fmri_output_folder)



print("Processing fMRI data...")
for file in os.listdir(fmri_folder):
    if not file.endswith("xlsx"):
        continue
    print(f"    Processing file: {file}")
    fmri_file = os.path.join(fmri_folder, file)
    fmri_data = pd.read_excel(fmri_file, index_col=0)
    # Input is already Fisher Z-transformed and globally normalized

    fmri_array = fmri_data.to_numpy().copy()
    print("    Computing graph metrics...")
    # compute the rest of the graph metrics in a function and return a dataframe with them

    GraphMetrics_df = compute_graph_metrics(fmri_data)
    print("    Computing modularity...")

    # Search for gamma values that yield 5-7 communities (biologically plausible range)
    # Uses three-stage verification pipeline (see function docstring for details)
    # Returns stable community assignments via consensus clustering


    mod_dataframe= search_best_gamma_seed_combination(data_array=fmri_array, gamma_step = 0.0125)

    mod_dataframe.index = fmri_data.index
    output_file = "fMRI_" + file.split(".xlsx")[0] + "_GraphMetrics.xlsx"
    with pd.ExcelWriter(os.path.join(fmri_output_folder, output_file)) as writer:
        GraphMetrics_df.to_excel(writer, sheet_name="GraphMetrics")
        mod_dataframe.to_excel(writer, sheet_name="Modularity")




print("Processing fNIRS data...")
# for fNIRS
for file in os.listdir(fnirs_folder):
    if not file.endswith("xlsx"):
        continue
    print(f"    Processing file: {file}")
    fnirs_file = os.path.join(fnirs_folder, file)
    sheet_names = pd.ExcelFile(fnirs_file).sheet_names
    output_file = "fNIRS_" + file.split(".xlsx")[0] + "_GraphMetrics.xlsx"
    with pd.ExcelWriter(os.path.join(fnirs_output_folder, output_file)) as writer:
        for sheet in sheet_names:
            fnirs_data = pd.read_excel(fnirs_file, index_col=0, sheet_name=sheet)
            # Input is already Fisher Z-transformed and globally normalized
            fnirs_array = fnirs_data.to_numpy().copy()
            print(f"    Processing sheet: {sheet}")
            print("    Computing graph metrics...")
            GraphMetrics_df = compute_graph_metrics(fnirs_data)
            print("    Computing modularity...")


            mod_dataframe= search_best_gamma_seed_combination(data_array=fnirs_array, gamma_step = 0.0125)


            mod_dataframe.index = fnirs_data.index

            # Write the DataFrame to the Excel file
            GraphMetrics_df.to_excel(writer, sheet_name=f"{sheet}_GraphMetrics")
            mod_dataframe.to_excel(writer, sheet_name=f"{sheet}_Modularity")
