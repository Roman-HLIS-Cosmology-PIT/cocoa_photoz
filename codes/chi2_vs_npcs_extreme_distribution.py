import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import glob
import re

cocoa_path = os.getcwd()
assert cocoa_path.strip().endswith("Cocoa"), "Please, run this script from ./whatever/cocoa/Cocoa"
parser = argparse.ArgumentParser()
parser.add_argument('--model', type=str, required=True, help='Simulation name or tag')
parser.add_argument('--fiducial', type=str, required=True, help='Simulation name or tag')

args = parser.parse_args()
mod  = args.model
fid  = args.fiducial

rename={"sc1bd4":"DRM-D1",
        "sc1bd5":"DRM-D2",
        "sc1bd6":"DRM-D3",
        "sc1bd7":"DRM-D4"
        }

assert mod in rename.keys() and fid in rename.keys(), f"""The Model {mod} or the Fiducial {fid} were not found in the list {rename.keys()}."""

name = f"chi2_vs_npcs_Model{mod}_Fiducial{fid}_compute"

########################################################################
folder = f"{cocoa_path}/chi2_vs_npcs_weighted"  
pattern = os.path.join(folder, f"{name}_*.txt") 
all_files = glob.glob(pattern)
integer_files = [
    f for f in all_files
    if re.search(rf"{re.escape(name)}_(\d+)\.txt$", os.path.basename(f))
]
print(f"Found {len(integer_files)} files with integer suffixes")
########################################################################

input_filename = f"{cocoa_path}/chi2_vs_npcs_weighted/{name}"
output_filename = f"{cocoa_path}/chi2_vs_npcs_weighted/{name}_all.txt"
output_filename_unique = f"{cocoa_path}/chi2_vs_npcs_weighted/{name}_all_unique.txt"

open(output_filename,"w").close()

for i in range(1,len(integer_files)+1):
    with open(f"{input_filename}_{i}.txt","r") as f:
        lines = f.readlines()
        with open(output_filename,"a") as g:
            for line in lines:
                g.write(line)

data = np.genfromtxt(output_filename, delimiter=",")
_, unique_indices = np.unique(data[:, 0], return_index=True)
unique_data = data[sorted(unique_indices)]
np.savetxt(f"{output_filename_unique}", unique_data, delimiter=",", fmt=["%d", "%.10f"])

data = np.genfromtxt(f"{output_filename_unique}", delimiter=",", skip_header=1)
chi2_values = data[:, 1]

##### find extreme #####
chi2_max = np.max(chi2_values)
print(f"The maximum chi2 found so far is {chi2_max}")

# Plot
def plot():
    plt.figure(figsize=(8, 5))
    plt.hist(chi2_values, bins='auto', alpha=0.7, edgecolor='black')
    # plt.yscale("log")
    plt.title(f"Model: {rename[mod]}, Fiducial: {rename[fid]}", fontsize=14)
    plt.xlabel("χ²", fontsize=12)
    plt.ylabel("Frequency of extreme scenarios", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.4)

    plt.tight_layout()
    plt.savefig("test.pdf")
    return None