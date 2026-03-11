# Light on Broken Networks: Resting-State fNIRS as a Tool for Connectivity Mapping

**Author:** Foivos Kotsogiannis
**Email:** f.kotsogiannisteftsoglou@alumni.maastrichtuniversity.nl
**Last Updated:** January 9, 2026
**DOI** https://doi.org/10.64898/2026.03.06.710143

This repository contains code for:

Foivos Kotsogiannis et al. (2026). Light on Broken Networks: Resting-State fNIRS as a Tool for Connectivity Mapping

bioRxiv. https://doi.org/10.64898/2026.03.06.710143

**Please cite the original paper on bioRxiv**


---

## Overview

This project investigates resting-state functional connectivity (RSFC) in fNIRS data by comparing it with gold-standard fMRI connectivity patterns from the Human Connectome Project (HCP). The pipeline processes both fNIRS (from the MULPA dataset) and fMRI (from HCP) data, computes correlation matrices using bivariate and partial correlation methods, performs graph theoretical analyses, and evaluates the spatial correspondence between modalities.

## Prerequisites

### Software Requirements
- **Python 3.x** with packages: MNE-Python, NumPy, Pandas, SciPy, Matplotlib, Seaborn, PySide6, NetworkX
- **FSL (FMRIB Software Library)** for fMRI preprocessing and ROI extraction
- **Satori** (Brain Innovation) for fNIRS preprocessing
- **BrainNet Viewer** (optional) for 3D network visualization

### Data Requirements
- **fNIRS Data:** Preprocessed MULPA dataset in .snirf format (resting-state)
- **fMRI Data:** HCP resting-state fMRI data (3T, preprocessed)
- **Template Files:** MNI152 brain template, channel coordinates, ROI definitions

---

## Pipeline Overview

The complete pipeline consists of 6 major stages:

1. **fMRI Data Handling** - Download and organize HCP data
2. **Channel-to-ROI Mapping** - Create fMRI ROIs matching fNIRS channel locations
3. **Correlation Matrix Computation** - Calculate connectivity matrices for both modalities
4. **Statistical Comparison** - Compare fNIRS and fMRI connectivity patterns
5. **Graph Metrics Analysis** - Compute network properties
6. **Network Visualization** - Prepare data for BrainNet Viewer

---

## 1. fMRI Data Preparation

**Purpose:** Download, organize, and prepare HCP resting-state fMRI data for analysis

**Script Folder:** `1_fMRI_DataPreparation/1_HandlingHCPDownload/`

### 1.1 Download HCP Data
- Download HCP resting-state fMRI data (3T) from the HCP database
- Target file: `rfMRI_REST*_*_hp2000_clean_rclean_tclean.nii.gz` (preprocessed, cleaned data)
- Original structure: `SubjID/MNINonLinear/Results/rfMRI_REST*/`

### 1.2 Organize Folder Structure

**Execute:** `1_HCP_FolderClean.py`

**Function:** Reorganizes downloaded HCP data into simplified folder structure

**Input:**
- HCP folder with complex nested structure
- Multiple REST runs per subject (REST1_LR, REST1_RL, REST2_LR, REST2_RL)

**Output:**
- Simplified structure: `Subjects/SubjID_Rest3T.nii.gz`
- Single combined or selected run per subject

**Note:** Target: NIfTI file for LP run from session 1

## 2. Channel-to-ROI Mapping

**Purpose:** Create fMRI ROIs (spheres) that match fNIRS channel locations for direct comparison

**Script Folder:** `1_fMRI_DataPreparation/2_ChannelToROI/`

This section creates spherical ROIs in MNI space centered at fNIRS channel locations, then extracts fMRI time series from these ROIs.

### 2.1 Generate ROI Spheres

**Script Folder:** `1_FSL_SphereCreation/`

**Execute (in order):**
1. `1_sphere_creation.sh` - Creates individual spheres for each fNIRS channel
2. `Combine_SingSpheres.sh` - Combines individual spheres into single mask

**Method:**
- Uses FSL's `fslmaths` to create spherical ROIs
- Sphere radius: 5mm
- Coordinates: MNI152 space matching fNIRS montage

**Input:**
- fNIRS channel coordinates in MNI space (from fOLD toolbox)
- MNI152 template brain

**Output:**
- Individual sphere masks: `ROI#_Channnel#_point.nii.gz`
- Combined mask file with all ROIs

**Requirements:**
- FSL / bash


### 2.2 Extract ROI Time Series

**Script Folder:** `2_ExtractROI/`

**Execute (in order):**

1. **`1_CombineSpheres.py`**
   - Combines individual sphere masks into single multi-ROI mask
   - Assigns unique integer values to each ROI
   - **Output:** Combined ROI mask file

2. **`2_ExtractBoldFromROI.py`**
   - Extracts average BOLD time series from each ROI
   - Processes all subjects
   - **Input:** 
     - fMRI data files (.nii.gz)
     - Combined ROI mask
   - **Output:** 
     - Excel files: `SubjID_ROI_Timeseries.xlsx`
     - Format: Columns = ROIs, Rows = Time points

**Processing Details:**
- Computes mean BOLD signal across all voxels within each spherical ROI
- Time series length depends on HCP scan duration (~15 minutes)
- Saves as pandas DataFrame for easy correlation analysis

---

## 3. Correlation Matrix Computation

**Purpose:** Calculate functional connectivity matrices using bivariate and partial correlation methods for both fNIRS and fMRI data

### 3.1 fMRI Correlation Matrices

**Script Folder:** `2_fMRI_CorrelationMatrices/`

**Execute (in order):**

1. **`1_StandardCor.py`** - Standard (Bivariate) Correlation
   - Computes Pearson correlation between all ROI pairs
   - **Input:** ROI time series Excel files
   - **Output:** Subject-level correlation matrices (ROI × ROI)
   - File suffix: `_StandardCorrelation.xlsx`

2. **`2_PartialCor.py`** - Partial Correlation
   - Computes partial correlation controlling for all other ROIs
   - Removes indirect connectivity influences
   - **Input:** ROI time series Excel files
   - **Output:** Subject-level partial correlation matrices (ROI × ROI)
   - File suffix: `_PartialCorrelation.xlsx`

3. **`3_AverageCorM.py`** - Group-Level Averaging
   - Applies Fisher Z-transformation to individual correlation matrices
   - Averages across subjects

   - **Input:** Individual correlation matrices (fMRI + fNIRS)
   - **Output:** 
     - Group-averaged correlation matrix
     - File: `fMRI_Group_Average_[StandardCor/PartialCor].xlsx`

4. **`4_PlotHeatAver.py`** - Visualization
   - Creates heatmap visualizations of group-average matrices
   - **Input:** Group-averaged correlation matrices (fMRI + fNIRS)
   - **Output:** Publication-quality heatmap figures (.png)

**STOP HERE:** Before proceeding to Histograms, compute fNIRS correlation matrices (see Section 3.2 below). Complete both bivariate and partial correlation for fNIRS data, then return here to continue.

5. **`5_CreateHisto.py`** - Distribution Analysis
   - Creates histograms of correlation values
   - Analyzes distribution properties
   - **Input:** Group-averaged matrices (from both fMRI/fNIRS Average Correlation Matrices)
   - **Output:** Histogram plots

### 3.2 fNIRS Correlation Matrices

**Script Folders:** `3_fNIRS_CorrelationMatrices/Biv Cor/` and `3_fNIRS_CorrelationMatrices/Partial Cor/`

**Method:** Similar pipeline to fMRI with adaptations for fNIRS channel data

#### Bivariate Correlation (`fNIRS/Biv Cor/`)

**Execute (in order):**

1. **`1_subject_level_cor.py`**
   - Loads preprocessed .snirf files
   - Removes short-distance channels (detector > 28)
   - Computes Pearson correlation between all channel pairs
   - **Input:** Preprocessed fNIRS data (.snirf format)
   - **Output:** Subject-level correlation matrices
   - File suffix: `_BivariateCor.xlsx`

2. **`2_AverageCorM.py`**
   - Fisher Z-transformation and averaging
   - **Output:** Group-average matrix

**Utility:** `utils.py` - Helper functions for fNIRS data handling

#### Partial Correlation (`fNIRS/Partial Cor/`)

**Execute (in order):**

1. **`1_subject_level_parcor.py`**
   - Computes partial correlation matrices
   - Controls for all other channels
   - **Input:** Preprocessed fNIRS data
   - **Output:** Subject-level partial correlation matrices

2. **`2_AverageCorM.py`**
   - Group-level averaging with Z-transformation

3. **`4_PlotHeatAver.py`**
   - Visualization of averaged matrices

---

## 4. Statistical Comparison of Correlation Matrices

**Purpose:** Compare fNIRS and fMRI connectivity patterns using edge-wise statistics and matrix-level metrics

**Script Folder:** `4_StatisticalComparison/`

### 4.1 Edge-wise Statistical Testing

**Script Folder:** `4_StatisticalComparison/1_EdgewiseTests/`

#### 4.1.1 Normality Testing

**Execute:**
- `1a_EdgewiseNormality_StandardCor.py` - Test normality for standard correlation
- `1c_EdgewiseNormality_PartialCor.py` - Test normality for partial correlation

**Function:**
- Tests each edge (connection) for normal distribution across subjects
- Uses Shapiro-Wilk test
- Determines appropriate statistical test (parametric vs. non-parametric)

**Output:**
- Normality test results per edge
- Summary statistics

#### 4.1.2 T-tests for Edge-wise Comparison

**Execute:**
- `1b_EdgewiseTtest_StandardCor.py` - T-tests for standard correlation edges
- `1d_EdgewiseTtest_PartialCor.py` - T-tests for partial correlation edges

**Function:**
- Compares each connection between fNIRS and fMRI
- independent t-tests because our fmri subjects are different from the fnirs subjects/ we have unpaired data
- Applies multiple comparison correction (FDR or Bonferroni)

**Output:**
- T-values and p-values for each edge
- Significant differences identified

#### 4.1.3 Visualization

**Execute:** `1e_Visualize_Edge_ttest.py`

**Function:**
- Creates visualizations of edge-wise comparison results
- Highlights significant differences
- Generates difference matrices (fNIRS - fMRI)

**Output:**
- Heatmaps showing statistical differences
- Color-coded significance levels

#### 4.1.4 TOST Equivalence Testing

**Execute:**
- `2a_TOST_Equivalence.py` - Two One-Sided Tests for equivalence
- `2b_Visualize_TOST.py` - Visualize TOST results

### 4.2 Matrix-Level Similarity

**Script Folder:** `4_StatisticalComparison/2_MatrixLevelTests/`

**Execute:** `2a_Mantel_CorrPermutation.py` - Correlation with Permutation Testing & Mantel Test
   
   **Purpose:** Compare fNIRS and fMRI group-averaged correlation matrices while accounting for non-independence of matrix elements
   
   **Methods Implemented:**
   
   a) **Pearson Correlation + Permutation Testing (Final Method)**
      - Computes Pearson correlation between upper triangle values of matrices
      - Performs permutation testing by shuffling matrix row/column labels
      - Preserves matrix dependency structure without distance conversion
      - **Advantage:** Direct correlation interpretation, no distance transformation needed
   
   b) **Mantel Test (Alternative/Verification)**
      - Converts correlation matrices to distance matrices (distance = 1 - correlation)
      - Computes correlation between distance matrices with permutation testing
      - Traditional method for comparing similarity matrices
   
   **Why Correlation + Permutation was chosen:**
   - Accounts for non-independence of matrix elements (like Mantel test)
   - Avoids conversion to distance matrices (preserves correlation interpretation)
   - More straightforward interpretation of similarity
   - Both methods yield equivalent statistical conclusions
   
   **Input:**
   - Group-averaged Fisher Z-transformed correlation matrices from fMRI and fNIRS
   
   **Output:**
   - `StandardCorrPermutation/` folder: Correlation coefficients and permutation-based p-values (HbO and HbR)
   - `MantelTest_Results.txt`: Mantel test results for verification
   - Permutations: 100,000 iterations

---
   - `StandardCorrPermutation/` folder: Correlation coefficients and permutation-based p-values (HbO and HbR)
   - `MantelTest_Results.txt`: Mantel test results for verification
   - Permutations: 100,000 iterations

**Process:**
- Computes Pearson correlation between upper triangle values of matrices
- Accounts for non-independence of matrix elements
- Verifies results with Mantel test

**Input:**
- Group-averaged Fisher Z-transformed correlation matrices from fMRI and fNIRS

**Output:**
- `StandardCorrPermutation/` folder: Correlation coefficients and permutation-based p-values (HbO and HbR)
- `MantelTest_Results.txt`: Mantel test results for verification

---

## 5. Graph Theoretical Analysis

**Purpose:** Characterize network topology using graph theory metrics

**Script Folder:** `5_GraphMetrics/`

### 5.1 Matrix Normalization

**Purpose:** Normalize correlation matrices to a common scale for cross-modal graph analysis comparisons

**Execute (in order):**

1. **`00a_NormalizeSubjctMatrices.py`** - Subject-Level Normalization
   
   **Process:**
   1. Apply Fisher Z-transformation to all subject correlation matrices (fMRI and fNIRS)
   2. Find global maximum absolute value across ALL subjects and modalities:
      - Search across all fMRI bivariate matrices -> find max
      - Search across all fNIRS bivariate matrices (HbO and HbR) -> find max
      - Determine global bivariate max from above
      - Repeat process for partial correlation matrices
   3. For each matrix:
      - Set diagonal values to the global max 
      - Divide entire matrix by the global max
      - Results in normalized values in range [0, 1]
   
   **Rationale:**
   - Fisher Z-space ensures proper handling of correlation distributions
   - Global normalization enables direct comparison between fMRI and fNIRS
   - Preserves relative strength of connections while standardizing scale
   
   **Input:** Subject-level correlation matrices (fMRI and fNIRS)
   
   **Output:** 
   - Normalized subject matrices in separate folders
   - Bivariate and partial handled separately with their own global maxima

2. **`00b_NormalizeGroupMatrices.py`** - Group-Level Normalization
   
   **Process:**
   1. Load group-averaged correlation matrices (already Fisher Z-transformed from averaging step)
   2. Find global maximum across group averages (fMRI, fNIRS HbO, fNIRS HbR)
   3. Normalize each group matrix by the global max
   
   **Note:** Does NOT re-apply Fisher Z-transformation (already applied during group averaging)
   
   **Input:** Group-averaged correlation matrices from Section 3
   
   **Output:** Normalized group matrices ready for graph analysis

**Important:**
- Bivariate and partial correlations use different global maxima
- Normalization preserves negative correlations (if present in Fisher Z-space)
- All subsequent graph metrics computed on these normalized matrices

### 5.2 Subject-Level Graph Metrics

**Execute:** `1_ComputeGraphMetrics_Subject.py`

**Function:**
- Computes node-level graph metrics for each subject's normalized connectivity matrix
- Operates on Fisher Z-transformed and globally normalized weighted graphs

**Metrics Computed (Node-Level Only):**

1. **Nodal Strength Metrics:**
   - Net Nodal Strength: Total weighted degree (sum of all connections)
   - Positive Nodal Strength: Sum of positive connections only
   - Negative Nodal Strength: Sum of negative connections only
   - Absolute Nodal Strength: Sum of positive + absolute value of negative
   - Normalized Nodal Strength: Based on Rubinov & Sporns (2011) formula accounting for positive/negative balance

2. **Local Assortativity:**
   - Positive Assortativity: Local assortativity using positive weights
   - Negative Assortativity: Local assortativity using negative weights
   - Computed using brainconn library (compatible with signed weights)

3. **Local Efficiency:**
   - Uses BCT (Brain Connectivity Toolbox) algorithm
   - Requires positive weights in range [0,1], therefore the absolute was used

4. **Nodal Path Length:**
   - Average shortest path length per node
   - Distance metric: 1 - |z| (treats positive and negative equally)

5. **Betweenness Centrality:**
   - Node importance based on shortest paths passing through each node
   - Uses same distance metric as path length: 1 - |z|

**Input:**
- Normalized subject-level correlation matrices (output from `00a_NormalizeSubjctMatrices.py`)
- Matrices in Fisher Z-space, normalized to [0, 1] range
- Both bivariate and partial correlation matrices processed

**Output:**
- Excel files with node-level metrics per subject
- fNIRS: Separate sheets for HbO and HbR
- fMRI: Single sheet per file
- File suffix: `_GraphMetrics.xlsx`
- Folder structure mirrors input (1_StandardCor/, 2_PartialCor/)

**Requirements:**
- BCT (Brain Connectivity Toolbox) Python implementation
- brainconn library for signed weight metrics


### 5.3 Group-Level Graph Metrics with Clustering

**Execute:** `2_ComputeGraphMetrics_Clustering.py`

**Function:**
- Computes nodal-level graph metrics for group-averaged connectivity matrices (same metrics as Section 5.2)
- Performs community detection using consensus clustering approach
- Searches for stable 5-7 community structures (biologically plausible based on known cortical systems)
- Implements three-stage computational efficiency pipeline

**Graph Metrics (Same as Subject-Level):**
- Nodal strength variants (net, positive, negative, absolute, normalized)
- Local assortativity (positive/negative)
- Local efficiency using |z|
- Nodal path length using distance 1-|z|
- Betweenness centrality using distance 1-|z|

**Community Detection (Consensus Clustering):**
- **Algorithm:** Louvain modularity optimization with consensus clustering
- **Target:** Find stable community structures with 5, 6, and 7 communities
- **Gamma Parameter:** Searched from 0.3 to ~2.0 (step: 0.0125 for final analysis)
  - Gamma < 1 -> Larger communities
  - Gamma > 1 -> Smaller communities

**Three-Stage Efficiency Pipeline:**
1. **Quick Screen (50 iterations):** Rapidly identify promising gamma values
2. **Verification (100 iterations):** Confirm gamma produces 5-7 communities
3. **Full Consensus (1000 iterations):** Generate final stable assignments via consensus clustering

**Consensus Clustering Process:**
- Runs Louvain 1000 times per gamma with random seeds
- Builds co-assignment matrix tracking node pair co-occurrences
- Applies consensus algorithm to extract most stable community structure
- Single consensus assignment represents all 1000 stochastic runs

**Output Selection:**
- For each community count (5, 6, 7): stores gamma with highest median modularity Q
- Q measures within-community edge density above chance
- Rejects solutions with singleton communities (isolated nodes)

**Output Files:**
- Excel files with two sheets:
  - **GraphMetrics:** Node-level metrics (same as subject-level)
  - **Modularity:** Three columns with consensus community assignments
    - Column format: "Q={median_Q}|G={gamma}|UC={count}"
    - Community labels: A, B, C, ... (letters for readability)

---

## 6. Comparison of Graph Metrics

**Purpose:** Compare network properties between fNIRS and fMRI at both group and subject levels

**Script Folder:** `6_GraphMetricsComparison/`

### 6.1 Group-Level Metric Comparison

**Execute:** `1_CompareGroupGraphMetrics.py`

**Function:**
- Correlates group-averaged patterns of graph metrics between fNIRS (HbO/HbR) and fMRI
- Tests relationship strength using multiple correlation methods
- Assesses normality of metric distributions

**Metrics Compared:**
- Net Nodal Strength
- Normalized Nodal Strength
- Positive Assortativity
- Local Efficiency
- Nodal Path Length (1-|z|)
- Betweenness Centrality (1-|z|)

**Statistical Tests:**
1. **Shapiro-Wilk Test:** Checks normality of metric distributions (fMRI and fNIRS separately)
2. **Pearson Correlation:** Parametric correlation coefficient
3. **Spearman Correlation:** Non-parametric rank correlation
4. **Kendall's Tau:** Alternative non-parametric correlation

**Process:**
- Loads group-level graph metrics from Section 5.3 outputs
- For each metric, computes all three correlation types
- Records correlation coefficients, p-values, and normality test results
- Captures warnings/errors during computation | like in the case of betweenness centrality

**Input:**
- Group-level graph metrics from Section 5.3 outputs
- `5_GraphMetrics/` (group-level outputs)
- Bivariate and partial correlation metrics processed separately

**Output:**
- Excel files: `GraphMetricCorrelation_between_fMRI_and_fNIRS_{type}_{chromophore}.xlsx`
- Single sheet with columns:
  - Graph Metric name
  - Pearson r & p-value
  - Spearman rho & p-value
  - Kendall tau & p-value
  - Shapiro-Wilk p-values (fMRI and fNIRS)
  - Status (ok/warning/error)
  - Notes (any warnings or errors)

### 6.2 Subject-Level Metric Comparison

**Execute:** `2_CompareSubjectGraphMetrics.py`

**Function:**
- Node-by-node independent t-tests comparing fNIRS and fMRI across all subjects
- Computes effect sizes (Cohen's d) for each node and metric
- Applies FDR correction to control for multiple comparisons

**Statistical Method:**
- **Independent samples t-test** (two-tailed)
  - Rationale: Different subjects in fMRI (HCP) vs. fNIRS (MULPA) cohorts
  - Compares same node across different subject groups
- **Effect Size:** Cohen's d using pooled standard deviation
- **Multiple Comparison Correction:** FDR (Benjamini-Hochberg method)

**Process:**
1. Combine all subject-level metrics into dictionaries organized by metric type
2. For each metric and each node (N=82):
   - Extract values across all fMRI subjects
   - Extract values across all fNIRS subjects (HbO or HbR)
   - Compute independent t-test
   - Calculate Cohen's d: (mean_fMRI - mean_fNIRS) / pooled_SD
3. Apply FDR correction across all 82 nodes per metric
4. Count percentage of significant nodes after FDR correction

**Metrics Compared:**
- Net Nodal Strength
- Normalized Nodal Strength
- Positive Assortativity
- Local Efficiency
- Nodal Path Length
- Nodal Betweenness

**Input:**
- Subject-level graph metrics from Section 5.2 outputs
- Located in data output folders (not in script repository)
- Separate folders for bivariate (1_StandardCor) and partial (2_PartialCor)

**Output:**
- Excel files with separate sheets per metric:
  - `GraphMetric_nodal_ttest_fmri_vs_fnirs_bivariate_hbo.xlsx`
  - `GraphMetric_nodal_ttest_fmri_vs_fnirs_bivariate_hbr.xlsx`
  - `GraphMetric_nodal_ttest_fmri_vs_fnirs_partialcor_hbo.xlsx`
  - `GraphMetric_nodal_ttest_fmri_vs_fnirs_partialcor_hbr.xlsx`
- Each sheet contains (one row per node):
  - t-statistic
  - Raw p-value
  - FDR-corrected p-value
  - fMRI mean
  - fNIRS mean
  - Pooled standard deviation
  - Cohen's d effect size
  - Notes (e.g., warnings about zero variance)
- Summary text file: `ttest_summary.txt` with percentage of significant nodes per metric

---

## 7. Network Visualization and Module Analysis

**Purpose:** Prepare data for 3D brain network visualization, analyze modular structure overlap with canonical networks, and quantify cross-modal community correspondence

**Script Folder:** `7_NetworkVisualizationAnalysis/`

These scripts integrate data from both fNIRS and fMRI for visualization and cross-modal community analysis.

### 7.1 Modularity-Based Matrix Sorting

**Execute:** `1_MatrixSortedModularity.py`

**Function:**
- Reorders correlation matrix rows/columns based on detected community structure
- Groups nodes by module membership for visual clarity
- Creates publication-ready heatmaps with community boundaries highlighted

**Input:**
- Group-averaged correlation matrices (fMRI and fNIRS)
- Community assignments from Section 5.3 (Modularity sheet)

**Output:**
- Sorted correlation matrices (Excel files)
- Heatmap visualizations with module boundaries
- Color-coded by community membership

**Use Case:** Visualizing modular organization in connectivity matrices

### 7.2 BrainNet Viewer Data Preparation

**Execute (in order):**

1. **`2_BrainNetPrepDataFrame.py`**
   
   **Function:**
   - Creates master dataframes combining MNI coordinates, ROI mappings, community assignments, and nodal strength
   - Ensures compatibility between channel names across all data sources
   
   **Output:**
   - Excel files with complete dataframes for fMRI (bivariate/partial) and fNIRS (HbO/HbR for bivariate/partial)
   - Saved in `BrainNetFiles/1_RawDataFrames/` folders

2. **`3_BrainNetNodeFile.py`**
   
   **Function:**
   - Converts master dataframes to BrainNet Viewer `.node` file format
   - Maps community letters (A-G) to color codes (1-7)
   - Normalizes nodal strength to 2.5-5 range for node sizing
   
   **Output:**
   - `.node` files: x, y, z, color, size, label (tab-separated)
   - One file per gamma value with 5-7 communities
   - Saved in `BrainNetFiles/2_BrainNetNodes/` folders

### 7.3 Network Overlap Analysis with Canonical Atlases

**Purpose:** Quantify voxel overlap between detected communities and established Yeo functional networks

**Execute:**

1. **`4alt_7YEOCheckIfInMap.py`** - Yeo 7-Network Voxel Overlap
   
   **Function:**
   - Creates 5mm spheres around each node coordinate
   - Counts voxels within spheres that fall into each Yeo 7-network
   - Also performs exact-point checking (no sphere)
   - Aggregates counts per community to identify dominant networks
   
   **Yeo 7 Networks:**
   - Visual, Somatosensory, Dorsal Attention, Ventral Attention, Limbic, Frontoparietal, Default
   
   **Output:**
   - Two Excel files per dataset: `*_Sphere_Check_*.xlsx` and `*_Point_Check_*.xlsx`
   - Columns: Module, Network_Name, Count (voxels), Percentage of hits
   - Separate analyses for 5, 6, and 7 community solutions

2. **`4alt_17YEOCheckIfInMap.py`** - Yeo 17-Network Voxel Overlap
   
   **Function:**
   - Same methodology as 7-network version
   - Uses finer-grained 17-network parcellation
   
   **Output:**
   - Sphere and point-based voxel counts for 17 networks
   - More detailed network subdivisions (e.g., Default A/B/C/D)

**Important** For the final analysis the 17-canonical Yeo map was used


3. **`5_JaccardSubnetworkOverlap.py`** - Cross-Modal Community Correspondence
   
   **Function:**
   - Computes Jaccard similarity between all fMRI and fNIRS community pairs
   - Identifies best-matching community pairs across modalities
   - Performs permutation testing to assess statistical significance of overlaps
   - Analyzes 7-community solutions only (UC=7)
   
   **Method:**
   - **Jaccard Index:** J(A,B) = |A ∩ B| / |A ∪ B| (set-based node overlap)
   - Computed for all fMRI x fNIRS community pairs
   - For each fMRI community: selects fNIRS community with maximum Jaccard index
   - **Total Overlap Score:** Mean of maximum Jaccard indices across all fMRI communities
   
   **Permutation Test (100,000 iterations):**
   - Null hypothesis: Community overlap is no greater than chance
   - Procedure:
     - Fix fMRI community assignments
     - Randomly shuffle fNIRS node assignments while preserving community sizes
     - Compute max Jaccard for each fMRI community under null
     - Build null distribution (100K shuffles)
   - **Output:** One-sided p-values testing if observed overlap > chance level
   
   **Input:**
   - Community assignments from Section 5.3 (both modalities, 7-community solutions)
   - Processes bivariate and partial correlations separately
   
   **Output (Excel files):**
   - Columns: fMRI Community Label, fNIRS Community Label (best match), Jaccard Index, p_greater
   - One row per fMRI community + Total Overlap Score row
   - Separate files for HbO and HbR
   - Community labels remapped: A->1, B->2, ..., G->7


**Network Atlases Used:**
- Yeo 17-Network parcellation (Thomas Yeo et al., 2011)
- MNI152 template space

### 7.4 Publication Visualization

**Execute:**

**`99_BrainNetAltFig.py`** - Selective Node Visualization for Publication

**Function:**
- Generates BrainNet Viewer input files with selective node display
- Shows only fNIRS nodes that overlap with their respective fMRI communities
- Creates cleaner visualizations by hiding non-overlapping nodes
- **Primary script used for final paper figures**

---

## Important Notes

### Data Requirements:
- **fNIRS data must be preprocessed** 
- **fMRI data** should be HCP-preprocessed or equivalently cleaned data
- **Channel coordinates** in MNI space required for ROI creation
see file (`Channels to Brain Areas using fOLD_updated_droppedmissing.xlsx`) or the uncleaned version (Last Sheet:`Channels to Brain Areas using fOLD.xlsx`)

## Contact

For questions or issues with this pipeline:
- **Author:** Foivos Kotsogiannis
- **Email:** f.kotsogiannisteftsoglou@alumni.maastrichtuniversity.nl

---

## Acknowledgments

- Human Connectome Project (HCP) for fMRI data
- MULPA dataset collaborators for fNIRS data
- Brain Innovation (Satori software)


