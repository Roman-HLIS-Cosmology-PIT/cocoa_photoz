# Diogo Souza - Jul 14 2025
import glob
import numpy as np
import matplotlib.pyplot as plt
import h5py
from cobaya.yaml import yaml_load
from yaml import safe_load, dump
from cobaya.model import get_model
import sys, os
sys.path.insert(0, os.environ['ROOTDIR']+'/external_modules/code/CAMB/build/lib.linux-x86_64-'+os.environ['PYTHON_VERSION'])
import camb
from camb import model
import argparse
import contextlib
import shutil

@contextlib.contextmanager
def suppress_stdout_stderr():
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr

parser = argparse.ArgumentParser()
parser.add_argument('--model', type=str, required=True, help='Simulation name or tag')
parser.add_argument('--fiducial', type=str, required=True, help='Simulation name or tag')

args        = parser.parse_args()
mod         = args.model
fid         = args.fiducial
cocoa_path  = '/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa' # Change for your path 
start_point = 0
output      = f"{cocoa_path}/chi2_vs_npcs_weighted/chi2_vs_npcs_Model{mod}_Fiducial{fid}_compute_1.txt"

if os.path.exists(output):
    if os.path.getsize(output) > 0:
        # stop immediately
        print(f"❌ File already exists and is not empty: {output}")
        sys.exit(1)
    else:
        # do not exit — let program proceed
        print(f"🚨 File exists but is empty, continuing: {output}")

camb_path   = os.path.dirname(camb.__file__)
path_model  = f"{cocoa_path}/projects/roman_real/data/{mod}"  
sims_path   = f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{fid}/nz_samples_LHC0_pointZ_1e6_Roman_{fid}.h5'
nz_h5py     = h5py.File(sims_path,'r')
tot_sims    = int(nz_h5py[f'bin0'].shape[0])
Nz          = int(nz_h5py[f'bin0'].shape[1])
Nt          = 9
z           = np.array(nz_h5py['zbinsc'])

print(f'# Using CAMB {camb.__version__} installed at {camb_path}')
print(f"# {nz_h5py.keys()}")   # ['bin0', 'bin1', 'bin2', 'bin3', 'bin4', 'bin5', 'bin6', 'bin7', 'bin8', 'zbins', 'zbinsc']
print(f"# {nz_h5py[f'bin0']}") # shape (1000000, 46)
print(f"# Nt: {Nt}, Nz: {Nz}, tot_sims: {tot_sims}")
print(f"# MODEL SCENARION {mod}. i.e., the mean nz.")
print(f"# RESULTS ARE SAVED @ {output}")
print(f"# FIDUCIAL SCENARIO {fid}. i.e., will vary nzs one the million realizations.")
print()
print(f"# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PURPOSE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
print(f""" PURPOSE: Compute the χ² for the fixed mean nz from scenario {mod}   
           for EACH of the 1M realizations from scenario {fid}""")        
print(f"# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PURPOSE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
print()
print("#############################################################")
print("##################### START COMPUTATION #####################")
print("#############################################################")
print()
print(f"# simulation index, χ²")
print()

##################################################################

def init_yaml():
    info     = None
    template_yaml = f'{cocoa_path}/cocoa_photoz/yamls/roman_pca/modelvector_generator_template.yaml'
    info_txt = f'{cocoa_path}/cocoa_photoz/yamls/roman_pca/modelvector_generator_Model{mod}_Fiducial{fid}.yaml' # This yaml file MUST exist apriori, so we copy from a template
    
    shutil.copy(template_yaml,info_txt) # Create a instance of the template for a specific model and fiducial
    
    with open(info_txt) as f: 
        info = yaml_load(f).copy()
        
    # Adjust the paths for theory (camb) and likelihood (roman) codes
    del info["theory"]["camb"]["path"]
    info["packages_path"] = camb_path
    info["likelihood"]["roman_real.roman_real_cosmic_shear"]["path"] = \
        f"./external_modules/data/roman_real/{mod}"
    info["likelihood"]["roman_real.roman_real_cosmic_shear"]["path"] =\
          f'{cocoa_path}/external_modules/data/roman_real/{mod}'
    return info
info = init_yaml()

def init_files(sim):

    df0_template  = f"template_no_nz.dataset"  
    df1_template  = f"template_no_modelvector.dataset"  
    df0_temporary = f"{fid}_r{sim}_fidRR.dataset"  
    df1_temporary = f"{mod}_nbar_{fid}r{sim}_fidRR.dataset"  
    nz_sim_file   = f"projects/roman_real/data/{mod}/nz_{fid}r{sim}.nz" 

    return df0_template  , \
           df1_template  , \
           df0_temporary , \
           df1_temporary , \
           nz_sim_file

def nz_sim(sim):
    _1, _2, _3 , _4 , nz_sim_file  = init_files(sim)
    nz = np.zeros((Nz,Nt))

    for i in range(Nt):
        nz[:,i]=nz_h5py[f'bin{i}'][sim]/np.trapz(y=nz_h5py[f'bin{i}'][sim],x=z)

    nz = np.column_stack((z,nz))
    np.savetxt(nz_sim_file,nz)
    return None

def create_temporary_data_file(sim,eval):
    df0_template, df1_template, df0_temporary , df1_temporary , _  = init_files(sim)
    df01_template  = [df0_template , df1_template]
    df01_temporary = [df0_temporary, df1_temporary]

    with open(f"{path_model}/{df01_template[eval]}","r") as g: 
        lines = g.readlines()

    open(f"{path_model}/{df01_temporary[eval]}","w").close()

    with open(f"{path_model}/{df01_temporary[eval]}","w") as f: 
        match eval:
            case 0:
                for line in lines:
                    if line.strip().startswith("nz_lens_file"):
                        new_line = f"nz_lens_file = nz_{fid}r{sim}.nz\n"
                        f.write(new_line) 
                    elif line.strip().startswith("nz_source_file"):
                        new_line = f"nz_source_file = nz_{fid}r{sim}.nz\n"
                        f.write(new_line)
                    else:
                        f.write(line)
            case 1: 
                for line in lines:
                    if line.strip().startswith("data_file"):
                        new_line = f"data_file = {fid}_r{sim}_fidRR.modelvector\n"
                        f.write(new_line)
                    else:
                        f.write(line)    
    return df01_temporary[eval]

def chi2_vs_npcs(sim,eval):
    df_temporary = create_temporary_data_file(sim,eval)
    info["likelihood"]["roman_real.roman_real_cosmic_shear"]["data_file"] = df_temporary
    
    match eval:
        case 0:
            info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector"] = True
            info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector_file"] =\
                  f'{path_model}/{fid}_r{sim}_fidRR.modelvector'
            with suppress_stdout_stderr():
                model = get_model(info)   
                varied_points = {} # The yaml do not contain any varied parameter! 
                model.loglikes(varied_points,as_dict=False,return_derived=False)[0]
        case 1:
            info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector"] = False
            info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector_file"] = ''

            with open(output,"a") as f:
                with suppress_stdout_stderr():
                    model = get_model(info)
                    varied_points = {} # The yaml do not contain any varied parameter! 
                    loglikes_i = model.loglikes(varied_points,as_dict=False,return_derived=False)[0]
                chi2_i = -2*loglikes_i
                print(f"{sim}, {chi2_i}")     # write in the SLURM .out 
                f.write(f"{sim}, {chi2_i}\n") # write in the .txt file           
    return None        

def clear_files_for_r(r):
    # Define patterns that match the files you want to remove
    _,_,_3,_4,_5 = init_files(sim=r)
    patterns = [f"{path_model}/{_3}",
                f"{path_model}/{_4}",
                f"{path_model}/{_5}"]

    for pattern in patterns:
        for filepath in glob.glob(pattern):
            try:
                os.remove(filepath)
                # print(f"Deleted {filepath}")  # Uncomment for debugging
            except FileNotFoundError:
                pass  # OK: file might not have been created

for r in range(start_point,tot_sims):
    try:
        nz_sim(sim=r)
        for i in [0, 1]:
            chi2_vs_npcs(sim=r, eval=i)
    finally:
        clear_files_for_r(r)