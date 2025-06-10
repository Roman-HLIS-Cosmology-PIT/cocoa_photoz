# survey  = 'lsst_y1'
# survey  = 'roman_real'
survey  = 'des_y3'

import sys, platform, os
os.environ['OMP_NUM_THREADS'] = '8'
import matplotlib
import math
from matplotlib import pyplot as plt
import numpy as np
import euclidemu2
import scipy

if survey == 'lsst_y1':
    import cosmolike_lsst_y1_interface as ci
elif survey == 'roman_real':
    import cosmolike_roman_real_interface as ci
elif survey == 'des_y3':
    import cosmolike_des_y3_interface as ci

print('------------------------------')
print('cosmolike interface used', ci)
print('------------------------------')

from getdist import IniFile
import itertools
import iminuit
import functools
print(sys.version)
print(os.getcwd())

# GENERAL PLOT OPTIONS
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['xtick.bottom'] = True
matplotlib.rcParams['xtick.top'] = False
matplotlib.rcParams['ytick.right'] = False
matplotlib.rcParams['axes.edgecolor'] = 'black'
matplotlib.rcParams['axes.linewidth'] = '1.0'
matplotlib.rcParams['axes.labelsize'] = 'medium'
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['grid.linewidth'] = '0.0'
matplotlib.rcParams['grid.alpha'] = '0.18'
matplotlib.rcParams['grid.color'] = 'lightgray'
matplotlib.rcParams['legend.labelspacing'] = 0.77
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.format'] = 'pdf'
matplotlib.rcParams['text.usetex'] = False

CAMBAccuracyBoost = 1.1
non_linear_emul = 2
CLprobe="xi"

path_cocoa="/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/"
titles = {'lsst_y1': 'LSST Y1 Source', 'roman_real': 'Roman Real Y1 Source', 'des_y3': 'DES Y3 Source'}

IA_model = 0
IA_redshift_evolution = 3

def survey_params(survey):
    if survey == 'lsst_y1':
        path= path_cocoa+f"external_modules/data/{survey}"
        data_file="lsst_y1_M1_GGL0.05.dataset"
        return path, data_file, IA_model, IA_redshift_evolution
    
    elif survey == 'roman_real':
        path= path_cocoa+f"external_modules/data/{survey}"
        data_file="example1.dataset"
        return path, data_file, IA_model, IA_redshift_evolution

    elif survey == 'des_y3':
        path= path_cocoa+f"external_modules/data/{survey}"
        data_file="des_y3_real.dataset"
        return path, data_file, IA_model, IA_redshift_evolution
    
    else:
        print("\n-----CRITICAL-----")
        print(f'\nSurvey "{survey}" not found\n')
        print("-----CRITICAL-----\n")

(path, data_file, IA_model, IA_redshift_evolution) = survey_params(survey=survey)

# NEED TO START COSMOLIKE 

ini = IniFile(os.path.normpath(os.path.join(path, data_file)))

print()
print("============")
print(path)
print(ini)
print("============")
print()

ci.initial_setup()
ci.init_accuracy_boost(1.0, 1.0, int(1))
ci.init_cosmo_runmode(is_linear = False)
if survey=='roman_real':
    ggl_exclude = [[6,0],[7,0],[7,1]]
    ggl_exclude = np.array(ggl_exclude).flatten()
    ci.init_ggl_exclude(ggl_exclude)
else: None    
ci.init_redshift_distributions_from_files(
      lens_multihisto_file=ini.relativeFileName('nz_lens_file'),
      lens_ntomo=int(ini.int("lens_ntomo")), 
      source_multihisto_file=ini.relativeFileName('nz_source_file'),
      source_ntomo=int(ini.int("source_ntomo")))


ci.init_IA( ia_model = int(IA_model), 
            ia_redshift_evolution = int(IA_redshift_evolution))

# Init Cosmolike
ci.init_probes(possible_probes = CLprobe)

ci.init_binning(int(ini.int("n_theta")), 
                ini.float("theta_min_arcmin"), 
                ini.float("theta_max_arcmin"))

ci.init_data_real(ini.relativeFileName('cov_file'), 
                  ini.relativeFileName('mask_file'), 
                  ini.relativeFileName('data_file'))

#~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~
def get_important_dimensions(twopt):
    s_ntomo = ini.int("source_ntomo")
    l_ntomo = ini.int("lens_ntomo")
    n_θ = ini.int("n_theta")
    n_θmin = ini.float("theta_min_arcmin")
    n_θmax = ini.float("theta_max_arcmin")
    ξp_dim  = int(( s_ntomo*(s_ntomo+1) / 2 ) * n_θ)
    ξm_dim  = int((s_ntomo*(s_ntomo+1) / 2 ) * n_θ)
    gammat_dim = int((s_ntomo * l_ntomo ) * n_θ)
    w_dim = int(s_ntomo * n_θ)
    print(\
        f"""
        DIMENSION INFO START\n
        Data vector dimension order defined at: https://github.com/CosmoLike/cocoa_roman_real/blob/3a5fd28ec7f0c60e093990177bf3f70d58649517/data/calculate_mask.py#L112C5-L112C58\n
        Data vector dimension order is: [ξp, ξm, γt, w]
        Source N tomo: {s_ntomo}
        Lens N tomo: {l_ntomo}
        N θ: {n_θ} ; θ ∈ [{n_θmin, n_θmax}] arcmin
        dim ξ+: {ξp_dim}
        dim ξ-: {ξm_dim}
        dim γt: {gammat_dim}
        dim w: {w_dim}
        Total DV dim: {ξp_dim+ξm_dim+gammat_dim+w_dim}
            Cross check: ci.get_dv_masked().shape: {np.array(ci.get_dv_masked()).shape},
        COV dim: ({ξp_dim+ξm_dim+gammat_dim+w_dim},{ξp_dim+ξm_dim+gammat_dim+w_dim}).
            Cross check: get_dv_masked().shape: {np.array(ci.get_cov_masked()).shape},
        Fisher dim: ({s_ntomo}xlen(nz),dv_dim) @ (dv_dim,dv_dim) @ (dv_dim,{s_ntomo}xlen(nz)) = ({s_ntomo}xlen(nz),{s_ntomo}xlen(nz))\n
        DIMENSION INFO END
        """
    )
    twopts = {"xip":ξp_dim,"xim":ξm_dim,"gammat":gammat_dim,"w":w_dim}
    return twopts[twopt]

def get_fisher_matrix(dtwoptdn_relative_path,twopt="xip"):
    twopt_dim = get_important_dimensions(twopt)
    print()
    print()
    print('TEST',ci.get_inv_cov_masked())
    print()
    print()
    inv_cov_masked = np.array(ci.get_inv_cov_masked())
    dxipdn = np.genfromtxt(path_cocoa+dtwoptdn_relative_path)
    inv_cov_xip = inv_cov_masked[0:twopt_dim,0:twopt_dim]
    fisher_mat = dxipdn @ inv_cov_xip @ dxipdn.T
    return fisher_mat

dtwoptdn_relative_path = f"cocoa_photoz/results_jacobian/{survey}/test_central_difference/test_central_difference_eps0.0001.txt"
# x = get_fisher_matrix(dtwoptdn_relative_path,twopt="xip")
# print(x.shape)

def plot_fisher_matrix(dtwoptdn_relative_path,twopt="xip"):

    fisher_mat = get_fisher_matrix(dtwoptdn_relative_path,twopt=twopt)

    im = plt.imshow(fisher_mat,cmap='seismic',vmin=-1,vmax=+1)
    plt.colorbar(im)
    plt.title(f'Fisher matrix: {titles[survey]}')
    plt.savefig(f'./plot_fisher_matrix__plot_fisher_matrix_{twopt}__{survey}.pdf')
    return 0

plot_fisher_matrix(dtwoptdn_relative_path,twopt="xip")
