# Diogo Souza - Jul 14 2025
import shutil
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
parser.add_argument('--realization', type=int, required=True, help='Simulation name or tag')

args        = parser.parse_args()
mod         = args.model
fid         = args.fiducial
realization = args.realization
camb_path   = os.path.dirname(camb.__file__)
cocoa_path  = '/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa' # Change for your path 
path_model  = f"{cocoa_path}/projects/roman_real/data/{mod}"  
sims_path   = f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{fid}/nz_samples_LHC0_pointZ_1e6_Roman_{fid}.h5'
nz_h5py     = h5py.File(sims_path,'r')
Nz          = int(nz_h5py[f'bin0'].shape[1])
Nt          = 9
z           = np.array(nz_h5py['zbinsc'])

print(f'# CHECK EXTREME Using CAMB {camb.__version__} installed at {camb_path}')
print(f"# CHECK EXTREME {nz_h5py.keys()}")   # ['bin0', 'bin1', 'bin2', 'bin3', 'bin4', 'bin5', 'bin6', 'bin7', 'bin8', 'zbins', 'zbinsc']
print(f"# CHECK EXTREME {nz_h5py[f'bin0']}") # shape (1000000, 46)
print(f"# CHECK EXTREME Nt: {Nt}, Nz: {Nz}")
print(f"# ATTENTION ATTENTION")
print(f"# TESTING REALIZATION: {realization} FOR MODEL {mod} AND FIDUCIAL {fid}")
print(f"# ATTENTION ATTENTION")
print("# CHECK EXTREME #############################################################")
print("# CHECK EXTREME ##################### START COMPUTATION #####################")
print("# CHECK EXTREME #############################################################")
print()

##################################################################

def nz_sim(sim):

    nz = np.zeros((Nz,Nt))

    for i in range(Nt):
        nz[:,i]=nz_h5py[f'bin{i}'][sim]/np.trapz(y=nz_h5py[f'bin{i}'][sim],x=z)

    nz = np.column_stack((z,nz))
    np.savetxt(f"projects/roman_real/data/{mod}/{fid}_r{realization}_check_extreme.nz",nz)
    print("\033[32mSTEP 1:\033[0m")
    print(f"NZ `{fid}_r{realization}_check_extreme.nz` SAVED AT `projects/roman_real/data/{mod}/`\n")
    return None

def create_files():
    df0_template   = f"template_no_nz_check_extreme.dataset"  
    df1_template   = f"template_no_modelvector_check_extreme.dataset"  
    df0_temporary  = f"{fid}_r{realization}_fidRR_check_extreme.dataset"  
    df1_temporary  = f"{mod}_nbar_{fid}_r{realization}_fidRR_check_extreme.dataset"  
    
    shutil.copy(f"{path_model}/{df0_template}",
                f"{path_model}/{df0_temporary}") # Create a instance of the template for a specific model and fiducial
    shutil.copy(f"{path_model}/{df1_template}",
                f"{path_model}/{df1_temporary}") # Create a instance of the template for a specific model and fiducial
    
    return df0_template,df1_template,df0_temporary,df1_temporary

df0_template, df1_template, df0_temporary, df1_temporary = create_files()

def create_temporary_data_file(sim,eval):
    df01_template  = [df0_template , df1_template , df1_template]
    df01_temporary = [df0_temporary, df1_temporary, df1_temporary]

    with open(f"{path_model}/{df01_template[eval]}","r") as g: 
        lines = g.readlines()

    open(f"{path_model}/{df01_temporary[eval]}","w").close()

    with open(f"{path_model}/{df01_temporary[eval]}","w") as f: 
        match eval:
            case 0:
                for line in lines:
                    if line.strip().startswith("nz_lens_file"):
                        new_line = f"nz_lens_file = {fid}_r{realization}_check_extreme.nz\n"
                        f.write(new_line) 
                    elif line.strip().startswith("nz_source_file"):
                        new_line = f"nz_source_file = {fid}_r{realization}_check_extreme.nz\n"
                        f.write(new_line)
                    else:
                        f.write(line)
            case 1:
                for line in lines:
                    if line.strip().startswith("data_file"):
                        new_line = f"data_file = {fid}_r{sim}_fidRR_check_extreme.modelvector\n"
                        f.write(new_line)
                    elif line.strip().startswith("nz_lens_file"):
                        new_line = f"nz_lens_file = {fid}_r{realization}_check_extreme.nz\n"
                        f.write(new_line) 
                    elif line.strip().startswith("nz_source_file"):
                        new_line = f"nz_source_file = {fid}_r{realization}_check_extreme.nz\n"
                        f.write(new_line)
                    else:
                        f.write(line)
            case 2:            
                for line in lines:
                    if line.strip().startswith("data_file"):
                        new_line = f"data_file = {fid}_r{sim}_fidRR_check_extreme.modelvector\n"
                        f.write(new_line)
                    else:
                        f.write(line)    
    return df01_temporary[eval]

def execute_cocoa_with_yaml(yaml_info):
    model = get_model(yaml_info)   
    varied_points = {} # The yaml do not contain any varied parameter! 
    loglikes_i = model.loglikes(varied_points,as_dict=False,return_derived=False)[0]
    return loglikes_i

def chi2_vs_npcs(sim,eval):
    ### INITIALIZE THE YAML FILE - START
    info     = None
    info_txt = f'{cocoa_path}/cocoa_photoz/yamls/roman_pca/modelvector_generator_check_extreme.yaml'
    with open(info_txt) as f: 
        info = yaml_load(f).copy()

    # Adjust the paths for theory (camb) and likelihood (roman) codes
    del info["theory"]["camb"]["path"]
    info["packages_path"] = camb_path
    info["likelihood"]["roman_real.roman_real_cosmic_shear"]["path"] =\
          f'{cocoa_path}/external_modules/data/roman_real/{mod}'
    ### INITIALIZE THE YAML FILE - END

    df_temporary = create_temporary_data_file(sim,eval)
    info["likelihood"]["roman_real.roman_real_cosmic_shear"]["data_file"] = df_temporary
    
    match eval:
        case 0:
            info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector"] = True
            info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector_file"] =\
                  f'{path_model}/{fid}_r{sim}_fidRR_check_extreme.modelvector'
            with suppress_stdout_stderr():
                loglikes_i = execute_cocoa_with_yaml(info)
            chi2_i = -2*loglikes_i    
            df_case0 = info["likelihood"]["roman_real.roman_real_cosmic_shear"]["data_file"]
            with open(f"{path_model}/{df_case0}") as f0:
                lines0 = f0.readlines()
                for line0 in lines0:
                    if line0.strip().startswith("data_file"):
                        line0_temp1 = line0
                    elif line0.strip().startswith("nz_lens_file"):
                        line0_temp2 = line0
            print("\033[32mSTEP 2:\033[0m")
            print(f'DATA FILE @ YAML: {df_case0}') 
            print(line0_temp1.rstrip("\n") + " \033[31m(ARBITRARY DV)\033[0m")            
            print(line0_temp2.rstrip("\n") + " \033[31m(CORRECT NZ)\033[0m") 
            print(f'print_datavector: {info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector"]}')               
            print(f'print_datavector_file: {info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector_file"]}')               
            print(f"SIM: {sim}, χ²: {chi2_i}\n")  

        case 1:
            info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector"] = False
            info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector_file"] = ''

            with suppress_stdout_stderr():
                loglikes_i = execute_cocoa_with_yaml(info)
            chi2_i = -2*loglikes_i
            df_case1 = info["likelihood"]["roman_real.roman_real_cosmic_shear"]["data_file"]
            with open(f"{path_model}/{df_case1}") as f1:
                lines1 = f1.readlines()
                for line1 in lines1:
                    if line1.strip().startswith("data_file"):
                        line1_temp1 = line1
                    elif line1.strip().startswith("nz_lens_file"):
                        line1_temp2 = line1
            print("\033[32mSTEP 3:\033[0m")
            print(f'DATA FILE @ YAML: {df_case1}') 
            print(line1_temp1.rstrip("\n") + " \033[31m(CORRECT DV)\033[0m")            
            print(line1_temp2.rstrip("\n") + " \033[31m(CORRECT NZ)\033[0m") 
            print(f'print_datavector: {info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector"]}')               
            print(f"SIM: {sim}, χ²: {chi2_i}\n") 

        case 2:
            info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector"] = False
            info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector_file"] = ''

            with suppress_stdout_stderr():
                loglikes_i = execute_cocoa_with_yaml(info)
            chi2_i = -2*loglikes_i
            df_case2 = info["likelihood"]["roman_real.roman_real_cosmic_shear"]["data_file"]
            with open(f"{path_model}/{df_case2}") as f2:
                lines2 = f2.readlines()
                for line2 in lines2:
                    if line2.strip().startswith("data_file"):
                        line2_temp1 = line2  
                    if line2.strip().startswith("nz_lens_file"):
                        line2_temp2 = line2  
            print("\033[32mSTEP 4:\033[0m")
            print(f'DATA FILE @ YAML: {df_case2}') 
            print(line2_temp1.rstrip("\n") + " \033[31m(CORRECT DV)\033[0m")            
            print(line2_temp2.rstrip("\n") + " \033[31m(MEAN NZ)\033[0m")        
            print(f'print_datavector: {info["likelihood"]["roman_real.roman_real_cosmic_shear"]["print_datavector"]}')               
            print(f"SIM: {sim}, χ²: {chi2_i}\n")     
    return None        


nz_sim(sim=realization)
for i in [0,1,2]:
    chi2_vs_npcs(sim=realization,eval=i)