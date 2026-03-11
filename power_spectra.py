

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from scipy.signal import welch

# --- Parameters ---
func_file = r"C:\Users\foivo\Desktop\Ubuntu Sharing\fMRI data\Subjects\102109_Rest3TRecommended\102109\MNINonLinear\Results\rfMRI_REST1_LR\rfMRI_REST1_LR_hp2000_clean_rclean_tclean.nii.gz"
TR = 0.72  # Repetition time in seconds
fs = 1.0 / TR  # Sampling frequency

# --- Load 4D fMRI data ---
img = nib.load(func_file)
data = img.get_fdata()

# --- Select one voxel (x, y, z) ---
voxel_coords = (50, 60, 40)  # Change as needed
voxel_ts = data[voxel_coords[0], voxel_coords[1], voxel_coords[2], :]

# --- Compute and plot power spectrum for voxel ---
f_voxel, psd_voxel = welch(voxel_ts, fs=fs, nperseg=256)
plt.figure()
plt.semilogy(f_voxel, psd_voxel)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power spectral density")
plt.title(f"Voxel power spectrum at {voxel_coords}")
plt.tight_layout()
plt.show()

# --- Compute and plot power spectrum for global mean ---
brain_ts = data.mean(axis=(0, 1, 2))
f_global, psd_global = welch(brain_ts, fs=fs, nperseg=256)
plt.figure()
plt.semilogy(f_global, psd_global)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power spectral density")
plt.title("Global mean power spectrum")
plt.tight_layout()
plt.show()
