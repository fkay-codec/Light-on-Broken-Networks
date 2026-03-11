#!/usr/bin/env bash

# Combine all ROI sphere files into a single mask using FSL's fslmaths.
# This script searches for files ending with '5mmSphere.nii.gz' in the input directory,
# and adds them together to create a combined ROI mask.

# we replicate this command from Andy:  
#   fslmaths roi1 -add roi2 -add roi3 outputfile

input_dir="/home/mybox/Desktop/Analysis/OutPuts"   # Adjust to your ROI sphere files directory
output_file="${input_dir}/Combined_ROIs.nii.gz"    # Output file name

# Find all 5mm sphere files
sphere_files=($(ls "$input_dir"/*5mmSphere.nii.gz))

# Check if any sphere files were found
if [ ${#sphere_files[@]} -eq 0 ]; then
    echo "No sphere files found in $input_dir"
    exit 1
fi

# print the count of sphere_files found; count should be 82
echo "Found ${#sphere_files[@]} sphere files."
# Check if the count is exactly 82
if [ ${#sphere_files[@]} -ne 82 ]; then
    echo "Error: Expected 82 sphere files, but found ${#sphere_files[@]}. Exiting."
    exit 1
fi

# Start building the fslmaths command or roi1
cmd="fslmaths \"${sphere_files[0]}\""

# Add each subsequent sphere file
for ((i=1; i<${#sphere_files[@]}; i++)); do
    cmd+=" -add \"${sphere_files[$i]}\""
done

cmd+=" \"$output_file\""

echo "Running command:"
echo $cmd
#exit 
# Execute the command
eval $cmd

echo "Combined ROI mask saved to $output_file"
