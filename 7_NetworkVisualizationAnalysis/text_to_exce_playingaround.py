import pandas as pd
import os
import re

# Load your text file
txt_path = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\7_GraphMetricsGroupLevel\BrainNetFiles\1_RawDataFrames\fMRI_Bivariate_Corr_Master_DF_Modules_in_RSFC_Check.txt"
txt_name = os.path.basename(txt_path)
with open(txt_path, 'r') as f:
    lines = f.readlines()

data = []
current_gamma = None
current_module = None

for line in lines:
    line = line.strip()
    
    # Detect Gamma
    gamma_match = re.match(r'For (G=[0-9.]+)', line)
    if gamma_match:
        current_gamma = gamma_match.group(1)
        print(f"Processing {current_gamma}")
        # quit()
        continue

    # Detect Module
    module_match = re.match(r'For Module ([A-Z]+):', line)
    if module_match:
        current_module = 'Module ' + module_match.group(1)
        row = [current_gamma, current_module]
        counts = []
        continue

    # Detect map points
    map_match = re.match(r'Map: (\S+), Points inside: (\d+)', line)
    if map_match:
        counts.append(int(map_match.group(2)))
        continue

    # End of module section
    if line.startswith("Highest count"):
        row.extend(counts)
        data.append(row)

# Column names
columns = [
    'Gamma', 'Module',
    'rsfmrinetwork_anterior_cingulate_precun',
    'rsfmrinetwork_auditory',
    'rsfmrinetwork_default',
    'rsfmrinetwork_IFG_middle_temporal',
    'rsfmrinetwork_left_executive',
    'rsfmrinetwork_left_right_exec_combined_network',
    'rsfmrinetwork_motor',
    'rsfmrinetwork_parietal_association_cortex',
    'rsfmrinetwork_posterior_default',
    'rsfmrinetwork_right_executive',
    'rsfmrinetwork_salience',
    'rsfmrinetwork_supplementary_motor',
    'rsfmrinetwork_visual'
]

df = pd.DataFrame(data, columns=columns)
import os
# Save to Excel
path_to_save = os.path.dirname(txt_path)
file_name=os.path.join(path_to_save, f'{txt_name}_counts.xlsx')
df.to_excel(file_name, index=False)
print("Excel file created: fMRI_module_counts.xlsx")
