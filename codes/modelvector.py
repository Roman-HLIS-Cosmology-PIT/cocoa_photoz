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
import yaml
from types import SimpleNamespace


surveys    = ['lsst_y1','des_y3','roman_real','roman_scenarios']
survey     = surveys[2]
path_c     = '/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/'
# path       = f'{path_c}/external_modules/data/{survey}/sc1bd4_g/'
path       = f'{path_c}/external_modules/data/{survey}/roman_scenarios/'

data_files = {'lsst_y1':'lsst_y1_M1_GGL0.05.dataset',
              'des_y3':'des_y3_real.dataset',
              #'roman_real': 'example1.dataset',
              #'roman_real':'roman_sc1bd4_g.dataset',
              'roman_real':'roman_sc1bd4.dataset'}

data_file = data_files[survey]

############ GENERAL PARAMETERS ############
non_linear_emul       = 2
IA_model              = 0
IA_redshift_evolution = 3
CAMBAccuracyBoost     = 1.1
AccuracyBoost         = 1.0
k_per_logint          = 20
CLAccuracyBoost       = 1.0
CLIntegrationAccuracy = 1
kmax                  = 10

############ LOAD SPEC PARAMETERS ############
with open('./cosmo_survey_parameters.yaml','r') as f: 
    config = yaml.safe_load(f)
cp = config['cosmo'] # cosmological parameters
sp = config[survey]  # survey parameters

############ COSMOLOGICAL PARMAETERS ############
As_1e9 = cp['As_1e9']
ns     = cp['ns']
H0     = cp['H0']
omegab = cp['omegab']
omegam = cp['omegam']
w0pwa  = cp['w0pwa']
w      = cp['w']
mnu    = cp['mnu']

############ SURVEY PARMAETERS ############
M                 = sp['M']
A1                = sp['A1']
A2                = sp['A2']
BTA               = sp['BTA']
shear_photoz_bias = sp['shear_photoz_bias']
lens_photoz_bias  = sp['lens_photoz_bias']
galaxy_bias_b1    = sp['galaxy_bias_b1']
galaxy_bias_b2    = sp['galaxy_bias_b2']
galaxy_bias_bmag  = sp['galaxy_bias_bmag']
bias_model        = sp['bias_model']
ntheta            = sp['ntheta'] 
theta_min_arcmin  = sp['theta_min_arcmin']
theta_max_arcmin  = sp['theta_max_arcmin']

############ DERIVED PARAMETERS ############
CLAccuracyBoost = CLAccuracyBoost * AccuracyBoost
CLSamplingBoost = CLAccuracyBoost * AccuracyBoost
CLIntegrationAccuracy = max(0, CLIntegrationAccuracy + 5*(AccuracyBoost-1.0))

############ INIT COSMOLIKE ############
module_name = f'cosmolike_{survey}_interface'
ci = importlib.import_module(module_name)

print('\n---------------------------')
print(f'cosmolike interface is {ci}')
print('---------------------------\n')

ini = IniFile(os.path.normpath(os.path.join(path, data_file)))

ci.initial_setup()
ci.init_accuracy_boost(1.0, CLAccuracyBoost, int(CLIntegrationAccuracy))
ci.init_cosmo_runmode(is_linear = False)

if survey == 'roman_real':
    # ggl_exclude = [[6,0],[7,0],[7,1]]
    ggl_exclude = []
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

############ DEF CAMB FUNCTION ############
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

############ CALL CAMB FUNCTION ############
(log10k_interp_2D, z_interp_2D, lnPL, lnPNL, G_growth, z_interp_1D, chi) = get_camb_cosmology(omegam=omegam, 
    omegab=omegab, H0=H0, ns=ns, As_1e9=As_1e9, w=w, w0pwa=w0pwa, AccuracyBoost=AccuracyBoost, kmax=kmax,
    k_per_logint=k_per_logint, CAMBAccuracyBoost=CAMBAccuracyBoost)

############ FEED COSMOLIKE ############
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
ci.init_bias(bias_model=bias_model)
ci.set_nuisance_shear_calib(M = M)
ci.set_nuisance_shear_photoz(bias = shear_photoz_bias)
ci.set_nuisance_clustering_photoz(bias = lens_photoz_bias)
ci.set_nuisance_bias(B1 = galaxy_bias_b1,
                    B2 = galaxy_bias_b2,
                    B_MAG = galaxy_bias_bmag)
ci.set_nuisance_ia(A1 = A1, A2 = A2, B_TA = BTA)

############ ACCESS 2PCF ############

n_tomo = int(ini.int("source_ntomo"))

xi_pm   = np.array(ci.xi_pm_tomo())
xi_p    = xi_pm[0]
xi_m    = xi_pm[1]
gamma_t = np.array(ci.w_gammat_tomo())
w_theta = np.array(ci.w_gg_tomo())


cs_tomo = [(i,j) for i in range(n_tomo) for j in range(n_tomo) if j>=i] # tomo bins combination for cosmic-shear => xi^{ij}(theta)
xi_p = np.hstack([xi_p[:,t[0],t[1]] for t in cs_tomo])
xi_m = np.hstack([xi_m[:,t[0],t[1]] for t in cs_tomo])
print('dim[xi+], dim[xi-] =',xi_p.shape,xi_m.shape)

ggl_tomo = [(i,j) for i in range(n_tomo) for j in range(n_tomo)] # tomo bins combination for galaxy-galaxy lensing => gammat^{ij}(theta)
gamma_t = np.hstack([gamma_t[:,t[0],t[1]] for t in ggl_tomo])
print('dim[gamma_t] =',gamma_t.shape)

gg_tomo = [(i,j) for i in range(n_tomo) for j in range(n_tomo) if j==i] # tomo bins combination for galaxy-galaxy => w^i(theta)
w_theta = np.hstack([w_theta[:,t[0],t[1]] for t in gg_tomo])
print('dim[w]',w_theta.shape)

print('total dim =',xi_p.shape[0] + xi_m.shape[0] + gamma_t.shape[0] + w_theta.shape[0])

mv = np.hstack((xi_p,xi_m,gamma_t,w_theta)) # model vector combined - CosmoLike order [xi+,xi-,gammat,w]
print('dim[mv] =',mv.shape)

mv = list(enumerate(mv))
np.savetxt('lcdm_sc1bd4.modelvector',mv,fmt='%d %e')

def plot_xipm():
    """plot xi_+ and xi_- for different Roman scenarios - see the paper"""
    fig,ax = plt.subplots(1,2,figsize=(7,3))
    xi_p0 = xi_pm[0][:,0,0]
    xi_m0 = xi_pm[1][:,0,0]
    thetas = np.arange(25., 250., 15.)
    ax[0].plot(thetas,thetas*xi_p0,color='#1b5f6f',lw=3,ls='-',marker='o',label='REF-D1')
    ax[1].plot(thetas,thetas*xi_m0,color='#1b5f6f',lw=3,ls='-',marker='o')
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[0].set_xlabel(r'$\theta$',fontsize=13)
    ax[1].set_xlabel(r'$\theta$',fontsize=13)
    ax[0].set_ylabel(r'$\theta\times\xi_+^{00}$',fontsize=13)
    ax[1].set_ylabel(r'$\theta\times\xi_-^{00}$',fontsize=13)
    ax[0].legend(loc='best')
    plt.tight_layout()
    # plt.savefig('test.pdf')
    return None