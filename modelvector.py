import os
os.environ['OMP_NUM_THREADS'] = '8'
import matplotlib.pyplot as plt
import numpy as np
import euclidemu2
from getdist import IniFile
import importlib
import camb
from camb import model
import scipy


surveys = ['lsst_y1','des_y3','roman_real']
data_files = {'lsst_y1':'lsst_y1_M1_GGL0.05.dataset',
              'des_y3':'des_y3_real.dataset',
              'roman_real':'example1.dataset'}
survey  = surveys[0]
path    = '../external_modules/data/' + survey
data_file = data_files[survey]

#############################
# CLprobe = "xi"
# CLprobe = "3x2pt"
# CLprobe = "gammat"
CAMBAccuracyBoost = 1.1
non_linear_emul = 2
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
w0pwa = -0.9
w = -0.9
AccuracyBoost=1.0
k_per_logint = 20
CLAccuracyBoost = 1.0
CLIntegrationAccuracy = 1
kmax = 10
CLAccuracyBoost = CLAccuracyBoost * AccuracyBoost
CLSamplingBoost = CLAccuracyBoost * AccuracyBoost
CLIntegrationAccuracy = max(0, CLIntegrationAccuracy + 5*(AccuracyBoost-1.0))
#############################

module_name = f'cosmolike_{survey}_interface'
ci = importlib.import_module(module_name)

print('\n---------------------------')
print(f'cosmolike interface is {ci}')
print('---------------------------\n')

ini = IniFile(os.path.normpath(os.path.join(path, data_file)))

ci.initial_setup()
ci.init_accuracy_boost(1.0, 1.0, int(1))
ci.init_cosmo_runmode(is_linear = False)

if survey == 'roman_real':
    ggl_exclude = [[6,0],[7,0],[7,1]]
    ggl_exclude = np.array(ggl_exclude).flatten()
    ci.init_ggl_exclude(ggl_exclude)
else: pass

ci.init_redshift_distributions_from_files(
    lens_multihisto_file=ini.relativeFileName('nz_lens_file'),
    lens_ntomo=int(ini.int("lens_ntomo")), 
    source_multihisto_file=ini.relativeFileName('nz_source_file'),
    source_ntomo=int(ini.int("source_ntomo")))

# ci.init_probes(possible_probes = CLprobe)

# ci.external_nz_modeling=int(1)
# nz=np.genfromtxt('../external_modules/data/lsst_y1/lsst_y1_source.nz')
# ci.set_source_sample(nz)

################################################
def get_camb_cosmology(omegam = omegam, omegab = omegab, H0 = H0, ns = ns, 
                       As_1e9 = As_1e9, w = w, w0pwa = w0pwa, AccuracyBoost = 1.0, 
                       kmax = 10, k_per_logint = 20, CAMBAccuracyBoost=1.1):

    As = lambda As_1e9: 1e-9 * As_1e9
    wa = lambda w0pwa, w: w0pwa - w
    omegabh2 = lambda omegab, H0: omegab*(H0/100)**2
    omegach2 = lambda omegam, omegab, mnu, H0: (omegam-omegab)*(H0/100)**2-(mnu*(3.046/3)**0.75)/94.0708
    omegamh2 = lambda omegam, H0: omegam*(H0/100)**2

    CAMBAccuracyBoost = CAMBAccuracyBoost*AccuracyBoost
    kmax = max(kmax/2.0, kmax*(1.0 + 3*(AccuracyBoost-1)))
    k_per_logint = max(k_per_logint/2.0, int(k_per_logint) + int(3*(AccuracyBoost-1)))
    extrap_kmax = max(max(2.5e2, 3*kmax), max(2.5e2, 3*kmax) * AccuracyBoost)

    z_interp_1D = np.concatenate( (np.concatenate( (np.linspace(0,2.0,1000),
                                                    np.linspace(2.0,10.1,200)),
                                                    axis=0
                                                 ),
                                   np.linspace(1080,2000,20)),
                                   axis=0)
    
    z_interp_2D = np.concatenate((np.linspace(0, 2.0, 95), np.linspace(2.25, 10, 5)),  axis=0)

    log10k_interp_2D = np.linspace(-4.2, 2.0, 1200)

    pars = camb.set_params(H0=H0, 
                           ombh2=omegabh2(omegab, H0), 
                           omch2=omegach2(omegam, omegab, mnu, H0), 
                           mnu=mnu, 
                           omk=0, 
                           tau=0.06,  
                           As=As(As_1e9), 
                           ns=ns, 
                           halofit_version='takahashi', 
                           lmax=10,
                           AccuracyBoost=CAMBAccuracyBoost,
                           lens_potential_accuracy=1.0,
                           num_massive_neutrinos=1,
                           nnu=3.046,
                           accurate_massive_neutrino_transfers=False,
                           k_per_logint=k_per_logint,
                           kmax = kmax);
    
    pars.set_dark_energy(w=w, wa=wa(w0pwa, w), dark_energy_model='ppf');    
    
    pars.NonLinear = model.NonLinear_both
    
    pars.set_matter_power(redshifts = z_interp_2D, kmax = kmax, silent = True);
    results = camb.get_results(pars)
    
    PKL  = results.get_matter_power_interpolator(var1="delta_tot", var2="delta_tot", nonlinear = False, 
                                                 extrap_kmax = extrap_kmax, hubble_units = False, k_hunit = False);
    
    PKNL = results.get_matter_power_interpolator(var1="delta_tot", var2="delta_tot",  nonlinear = True, 
                                                 extrap_kmax = extrap_kmax, hubble_units = False, k_hunit = False);
    
    lnPL = np.empty(len(log10k_interp_2D)*len(z_interp_2D))
    for i in range(len(z_interp_2D)):
        lnPL[i::len(z_interp_2D)] = np.log(PKL.P(z_interp_2D[i], np.power(10.0,log10k_interp_2D)))
    lnPL  += np.log(((H0/100.)**3)) 
    
    lnPNL  = np.empty(len(log10k_interp_2D)*len(z_interp_2D))
    if non_linear_emul == 1:
        params = { 'Omm'  : omegam, 
                   'As'   : As(As_1e9), 
                   'Omb'  : omegab,
                   'ns'   : ns, 
                   'h'    : H0/100., 
                   'mnu'  : mnu,  
                   'w'    : w, 
                   'wa'   : wa(w0pwa, w)
                 }
        kbt, bt = euclidemu2.get_boost( params, 
                                        z_interp_2D, 
                                        np.power(10.0, np.linspace( -2.0589, 0.973, len(log10k_interp_2D)))
                                      )
        log10k_interp_2D = log10k_interp_2D - np.log10(H0/100.)
        
        for i in range(len(z_interp_2D)):    
            lnbt = scipy.interpolate.interp1d(np.log10(kbt), np.log(bt[i]), kind = 'linear', 
                                              fill_value = 'extrapolate', 
                                              assume_sorted = True)(log10k_interp_2D)
            lnbt[np.power(10,log10k_interp_2D) < 8.73e-3] = 0.0
            lnPNL[i::len(z_interp_2D)]  = lnPL[i::len(z_interp_2D)] + lnbt
    elif non_linear_emul == 2:
        for i in range(len(z_interp_2D)):
            lnPNL[i::len(z_interp_2D)] = np.log(PKNL.P(z_interp_2D[i], np.power(10.0, log10k_interp_2D)))            
        log10k_interp_2D = log10k_interp_2D - np.log10(H0/100.)
        lnPNL += np.log(((H0/100.)**3))

    G_growth = np.sqrt(PKL.P(z_interp_2D,0.0005)/PKL.P(0,0.0005))
    G_growth = G_growth*(1 + z_interp_2D)
    G_growth = G_growth/G_growth[len(G_growth)-1]

    chi = results.comoving_radial_distance(z_interp_1D, tol=1e-4) * (H0/100.)

    return (log10k_interp_2D, z_interp_2D, lnPL, lnPNL, G_growth, z_interp_1D, chi)

(log10k_interp_2D, z_interp_2D, lnPL, lnPNL, G_growth, z_interp_1D, chi) = get_camb_cosmology(omegam=omegam, 
    omegab=omegab, H0=H0, ns=ns, As_1e9=As_1e9, w=w, w0pwa=w0pwa, AccuracyBoost=AccuracyBoost, kmax=kmax,
    k_per_logint=k_per_logint, CAMBAccuracyBoost=CAMBAccuracyBoost)

ci.init_accuracy_boost(1.0, CLAccuracyBoost, int(CLIntegrationAccuracy))

ci.init_binning(int(ntheta), theta_min_arcmin, theta_max_arcmin)

ci.set_cosmology(omegam = omegam, 
                    H0 = H0, 
                    log10k_2D = log10k_interp_2D, 
                    z_2D = z_interp_2D, 
                    lnP_linear = lnPL,
                    lnP_nonlinear = lnPNL,
                    G = G_growth,
                    z_1D = z_interp_1D,
                    chi = chi)
################################################

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
LSST_DZ_L1 = 0.00457604
LSST_DZ_L2 = 0.000309875
LSST_DZ_L3 = 0.00855907
LSST_DZ_L4 = -0.00316269
LSST_DZ_L5 = -0.0146753 
LSST_B1_1 = 1.72716
LSST_B1_2 = 1.65168
LSST_B1_3 = 1.61423
LSST_B1_4 = 1.92886
LSST_B1_5 = 2.11633

M = [LSST_M1, LSST_M2, LSST_M3, LSST_M4, LSST_M5]
shear_photoz_bias = [LSST_DZ_S1, LSST_DZ_S2, LSST_DZ_S3, LSST_DZ_S4, LSST_DZ_S5]
lens_photoz_bias = [LSST_DZ_L1, LSST_DZ_L2, LSST_DZ_L3, LSST_DZ_L4, LSST_DZ_L5]
galaxy_bias_b1 = [LSST_B1_1, LSST_B1_2, LSST_B1_3, LSST_B1_4, LSST_B1_5]
galaxy_bias_b2 = [0,0,0,0,0]
galaxy_bias_bmag = [0,0,0,0,0]
A1  = [LSST_A1_1, LSST_A1_2, 0, 0, 0]
A2  = [0, 0, 0, 0, 0]
BTA = [0, 0, 0, 0, 0]

ci.init_bias(bias_model=[0,0,0,1,0])
ci.set_nuisance_shear_photoz(bias = shear_photoz_bias)
ci.set_nuisance_clustering_photoz(bias = lens_photoz_bias)
ci.set_nuisance_bias(B1 = galaxy_bias_b1,
                    B2 = galaxy_bias_b2,
                    B_MAG = galaxy_bias_bmag)

################################################

xi_pm = np.array(ci.xi_pm_tomo())
gamma_t = np.array(ci.w_gammat_tomo())
w_theta = np.array(ci.w_gg_tomo())

# print(xi_pm.shape)
# print(w_gammat.shape)
# print(w_gg.shape)

print(gamma_t[0][0])