import sys, platform, os
os.environ['OMP_NUM_THREADS'] = '8'
import matplotlib
import math
from matplotlib import pyplot as plt
import numpy as np
# import euclidemu2
import scipy
import cosmolike_lsst_y1_interface as ci
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

path_cocoa="/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz"
path= "../../../external_modules/data/lsst_y1"
data_file="lsst_y1_M1_GGL0.05.dataset"

IA_model = 0
IA_redshift_evolution = 3

ntheta = 26 
theta_min_arcmin = 2.5 
theta_max_arcmin = 900

As_1e9 = 2.1
ns = 0.96605
H0 = 67.32
omegab = 0.04
omegam = 0.3
mnu = 0.06
LSST_DZ_S1 = 0.0414632
LSST_DZ_S2 = 0.00147332
LSST_DZ_S3 = 0.0237035
LSST_DZ_S4 = -0.0773436
LSST_DZ_S5 = -8.67127e-05
LSST_M1 = 0.0191832
LSST_M2 = -0.0431752
LSST_M3 = -0.034961
LSST_M4 = -0.0158096
LSST_M5 = -0.0158096
LSST_A1_1 = 0.606102
LSST_A1_2 = -1.51541
w0pwa = -0.9
w = -0.9

# NEED TO START COSMOLIKE 

ini = IniFile(os.path.normpath(os.path.join(path, data_file)))

ci.initial_setup()
ci.init_accuracy_boost(1.0, 1.0, int(1))
ci.init_cosmo_runmode(is_linear = False)

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
    inv_cov_masked = np.array(ci.get_inv_cov_masked())
    dxipdn = np.genfromtxt(path_cocoa+dtwoptdn_relative_path)
    inv_cov_xip = inv_cov_masked[0:twopt_dim,0:twopt_dim]
    fisher_mat = dxipdn @ inv_cov_xip @ dxipdn.T
    return fisher_mat

dtwoptdn_relative_path = "/results_jacobian/test_central_difference/test_central_difference_eps0.0001.txt"
# x = get_fisher_matrix(dtwoptdn_relative_path,twopt="xip")
# print(x.shape)

def plot_fisher_matrix(dtwoptdn_relative_path,twopt="xip"):

    fisher_mat = get_fisher_matrix(dtwoptdn_relative_path,twopt=twopt)

    im = plt.imshow(fisher_mat,cmap='seismic',vmin=-1,vmax=+1)
    plt.colorbar(im)
    plt.savefig('.figures/fisher_xip_lsst_y1.pdf')
    return 0

plot_fisher_matrix(dtwoptdn_relative_path,twopt="xip")
