"""
This script will create the same files as BrainNetParseData.py that are fed to BrainNet Viewer for vizualization.
The only difference is that this script the files will reflect our results from the paper.

Because of that this script is highly manual
"""

import numpy as np
import pandas as pd
import os
import shutil


out = r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\BrainNetAltFig"
if not os.path.exists(out):
    os.makedirs(out)

# copy the .node files you need to the out. We just keep the ones with 7 communities and also, it they exist in the out then we dont do it...
files2copy = [
    r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\7_GraphMetricsGroupLevel_Updated\BrainNetFiles\2_BrainNetNodes\fMRI_Bivariate_Q0p3961_G1p5125_UC7.node",
    r"c:\Users\foivo\Desktop\fNIRS_fMRI_Study\fMRI_data\7_GraphMetricsGroupLevel_Updated\BrainNetFiles\2_BrainNetNodes\fMRI_Partial_Q0p5215_G1p1125_UC7.node",
    r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\5_GraphMetricsGroupLevel_Updated\BrainNetFiles\2_BrainNetNodes\fNIRS_Bivariate_HbO_Q0p1650_G1p525_UC7.node",
    r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\5_GraphMetricsGroupLevel_Updated\BrainNetFiles\2_BrainNetNodes\fNIRS_Bivariate_HbR_Q0p1474_G1p2875_UC7.node",
    r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\5_GraphMetricsGroupLevel_Updated\BrainNetFiles\2_BrainNetNodes\fNIRS_Partial_HbO_Q0p3578_G1p1375_UC7.node",
    r"C:\Users\foivo\Desktop\fNIRS_fMRI_Study\fNIRS_data\5_GraphMetricsGroupLevel_Updated\BrainNetFiles\2_BrainNetNodes\fNIRS_Partial_HbR_Q0p3501_G1p1125_UC7.node",
]


for src in files2copy:
    base_name = os.path.basename(src)
    out_path = os.path.join(out, base_name)
    if os.path.exists(out_path):
        print(f"File {out_path} already exists, skipping...")
        continue
    if not os.path.exists(src):
        print(f"Source file {src} not found, skipping.")
        continue
    try:
        print(f"Copying {src} to {out_path}")
        shutil.copy2(src, out_path)
    except Exception as e:
        print(f"Failed to copy {src}: {e}")

# now things are more convenient, bcs they are grouped in one location

files_in_out = [f for f in os.listdir(out) if f.endswith('.node')]

#.node files are text with tab separated values with the following columns: X Y Z color size label
# we are going to take the color which is essentially the module label, and comment out lines that dont exist in the community. 
# so fMRI is fixed, and depending on the results we will comment out fNIRS nodes. 
# e.g., if fNIRS HbO comm1 is sig matched to fmri comm1, then we keep those nodes of fnirs hbo comm1 that are in the fmri comm1 and comment out those that are not there.





# lets make the variable to pass to the comment_out_nodes function
# these variables will be based on the results and are CONSTANT

bivariate_fmri_fnirs_hbo_comm_match = {
    1: 1,
    2: 2,
    3: 6,
    4: 5,
    5: 4,
    6: 7,
    7: 7,
}

bivariate_fmri_fnirs_hbr_comm_match = {
    1: 1,
    2: 2,
    3: 6,
    4: 6,
    5: 4,
    6: 7,
    7: 7,
}

partial_fmri_fnirs_hbo_comm_match = {
    1: 1,
    2: 2,
    3: 1,
    4: 6,
    5: 3,
    6: 5,
    7: 5,
}

partial_fmri_fnirs_hbr_comm_match = {
    1: 1,
    2: 2,
    3: 1,
    4: 4,
    5: 4,
    6: 7,
    7: 7,
}

def comment_out_nodes(fmri_node_file, fnirs_node_file, out_folder, fmri_fnirs_comm_match):
    """
    fmri_node_file: path to fmri .node file
    fnirs_node_file: path to fnirs .node file

    In this function, we tell which fmri-fnirs communities correspond to each other and comment out the fnirs that do not

    At the end we save the new fnirs .node file with a suffix _CommentedOut.node
    
    When we have 2 communities that correspond to the same fmri community, then we keep all nodes from both fnirs communities that match to that single fnirs community 
    """
    
    # Read fMRI nodes and build community -> labels mapping
    fmri_comm_labels = {}
    with open(fmri_node_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')


            if len(parts) >= 6:
                try:
                    community = int(float(parts[3]))  # color column
                    label = parts[5].strip()  # label column
                    if community not in fmri_comm_labels:
                        fmri_comm_labels[community] = set()
                    fmri_comm_labels[community].add(label)
                except:
                    continue
    
    # For each fNIRS community, determine which fMRI communities match to it
    # and build the set of allowed labels (union of all matched fMRI community labels)
    fnirs_allowed_labels = {}
    for fmri_comm, fnirs_comm in fmri_fnirs_comm_match.items():
        if fnirs_comm not in fnirs_allowed_labels:
            fnirs_allowed_labels[fnirs_comm] = set()

        # Add labels from this fMRI community
        if fmri_comm in fmri_comm_labels:
            # print(f"Matching fMRI community {fmri_comm} to fNIRS community {fnirs_comm} with labels: {fmri_comm_labels[fmri_comm]}\n")
            fnirs_allowed_labels[fnirs_comm].update(fmri_comm_labels[fmri_comm])
    
    # # Debug print
    # for comm, labels in fnirs_allowed_labels.items():
    #     print(f"fNIRS community {comm} allowed labels: {labels}\n")


    # Read fNIRS node file and process
    with open(fnirs_node_file, 'r') as f:
        lines = f.readlines()
    
    output_lines = []
    for line in lines:
        line = line.rstrip('\n\r')
        # Keep already commented or empty lines as is
        if line.strip().startswith('#') or not line.strip():
            output_lines.append(line)
            continue
        
        parts = line.split('\t')

        if len(parts) < 6:
            output_lines.append(line)
            continue
        
        try:
            community = int(float(parts[3]))  # color column
            label = parts[5].strip()  # label column
        except:
            output_lines.append(line)
            continue
        
        # Check if this fNIRS node should be kept
        if community in fnirs_allowed_labels:

            if label in fnirs_allowed_labels[community]:
                # Keep the node
                output_lines.append(line)
            else:
                # Comment out - label not in matched fMRI community
                output_lines.append('#' + line)
        else:
            # Community not matched at all, comment out
            output_lines.append('#' + line)
    
    # Create output filename
    base_name = os.path.basename(fnirs_node_file)
    name_without_ext = os.path.splitext(base_name)[0]
    output_filename = name_without_ext + '_CommentedOut.node'
    output_path = os.path.join(out_folder, output_filename)
    
    # Write the output file
    with open(output_path, 'w') as f:
        f.write('\n'.join(output_lines) + '\n')
    
    print(f"Created {output_path}")
    return output_path

fmri_bivariate_node = [f for f in files_in_out if f.startswith("fMRI_Bivariate") and f.endswith(".node")]
fnirs_bivariate_hbo_node = [f for f in files_in_out if f.startswith("fNIRS_Bivariate_HbO") and f.endswith(".node") and not f.endswith("_CommentedOut.node")]
fnirs_bivariate_hbr_node = [f for f in files_in_out if f.startswith("fNIRS_Bivariate_HbR") and f.endswith(".node") and not f.endswith("_CommentedOut.node")]

fmri_partial_node = [f for f in files_in_out if f.startswith("fMRI_Partial") and f.endswith(".node")]
fnirs_partial_hbo_node = [f for f in files_in_out if f.startswith("fNIRS_Partial_HbO") and f.endswith(".node") and not f.endswith("_CommentedOut.node")]
fnirs_partial_hbr_node = [f for f in files_in_out if f.startswith("fNIRS_Partial_HbR") and f.endswith(".node") and not f.endswith("_CommentedOut.node")]

comment_out_nodes(
    fmri_node_file=os.path.join(out, fmri_bivariate_node[0]),
    fnirs_node_file=os.path.join(out, fnirs_bivariate_hbo_node[0]),
    out_folder=out,
    fmri_fnirs_comm_match=bivariate_fmri_fnirs_hbo_comm_match
)
comment_out_nodes(
    fmri_node_file=os.path.join(out, fmri_bivariate_node[0]),
    fnirs_node_file=os.path.join(out, fnirs_bivariate_hbr_node[0]),
    out_folder=out,
    fmri_fnirs_comm_match=bivariate_fmri_fnirs_hbr_comm_match
)
comment_out_nodes(
    fmri_node_file=os.path.join(out, fmri_partial_node[0]),
    fnirs_node_file=os.path.join(out, fnirs_partial_hbo_node[0]),
    out_folder=out,
    fmri_fnirs_comm_match=partial_fmri_fnirs_hbo_comm_match
)
comment_out_nodes(
    fmri_node_file=os.path.join(out, fmri_partial_node[0]),
    fnirs_node_file=os.path.join(out, fnirs_partial_hbr_node[0]),
    out_folder=out,
    fmri_fnirs_comm_match=partial_fmri_fnirs_hbr_comm_match
)
