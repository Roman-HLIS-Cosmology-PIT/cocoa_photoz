# Example to illustrate how to use cocoa's new feature at commit/89d9300178dd8f645cdc4b1dd3e9e0b55845e760
# Example to illustrate how to use cocoa_lsst_y1's new feature at commit/1e7965bdce7ce568dcc8e6401c81bae201e0ba60
# Commit: "added external_nz_modeling to cosmolike generic interface"

import numpy as np
import cosmolike_lsst_y1_interface as ci
# import euclidemu2
import getdist
import os
import camb
import scipy
import itertools
import iminuit
import functools
import matplotlib.pyplot as plt

non_linear_emul = 2
CLprobe="xi"
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

path = "../../../external_modules/data/lsst_y1"
data_file = "lsst_y1_M1_GGL0.05.dataset"

def init_cosmolike(external_nz_modeling,mod_nz):
    ini = getdist.inifile.IniFile(
        os.path.normpath(os.path.join(path,data_file)))    

    ci.initial_setup()
    ci.external_nz_modeling = external_nz_modeling
    source_ntomo = ini.int("source_ntomo")
    lens_ntomo = ini.int("lens_ntomo")
    lens_file = ini.relativeFileName("nz_lens_file")

    (source_file, source_nz_local) = mod_nz()

    if (0 == external_nz_modeling):
        print(f"Old cosmolike feature: nz from file: {source_file}")
        ci.init_redshift_distributions_from_files(
            lens_multihisto_file=lens_file,
            lens_ntomo=int(lens_ntomo), 
            source_multihisto_file=source_file,
            source_ntomo=int(source_ntomo))

    elif (1 == external_nz_modeling):
        print("New cosmolike feature: nz from array: source_nz_local")
        ci.init_source_sample_size(int(source_ntomo))
        ci.init_lens_sample_size(int(lens_ntomo))
        ci.set_source_sample(source_nz_local)
        ci.init_ntomo_powerspectra()

    ci.init_accuracy_boost(1.0, 1.0, int(1))
    ci.init_cosmo_runmode(is_linear = False)
    ci.init_IA( ia_model = int(IA_model), 
            ia_redshift_evolution = int(IA_redshift_evolution))

    ci.init_probes(possible_probes = CLprobe)
    ci.init_binning(int(ntheta), theta_min_arcmin, theta_max_arcmin)
    ci.init_data_real(ini.relativeFileName('cov_file'), 
                  ini.relativeFileName('mask_file'), 
                  ini.relativeFileName('data_file'))

# def modified_nz():
#     spath = "../../external_modules/data/lsst_y1/lsst_y1_source"
#     x = np.genfromtxt(spath+".nz")
#     x[:,1] += .1
#     np.savetxt(spath+"_modified.nz",x)
#     return spath+"_modified.nz",x

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
    
    pars.NonLinear = camb.model.NonLinear_both
    
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

def xi(external_nz_modeling,mod_nz,ntheta = ntheta, 
    theta_min_arcmin = theta_min_arcmin, 
    theta_max_arcmin = theta_max_arcmin, 
    omegam = omegam, 
    omegab = omegab, 
    H0 = H0, 
    ns = ns, 
    As_1e9 = As_1e9, 
    w = w, 
    w0pwa = w0pwa,
    A1  = [LSST_A1_1, LSST_A1_2, 0, 0, 0], 
    A2  = [0, 0, 0, 0, 0],
    BTA = [0, 0, 0, 0, 0],
    shear_photoz_bias = [LSST_DZ_S1, LSST_DZ_S2, LSST_DZ_S3, LSST_DZ_S4, LSST_DZ_S5],
    M = [LSST_M1, LSST_M2, LSST_M3, LSST_M4, LSST_M5],
    baryon_sims = None,
    AccuracyBoost = 1.0, 
    kmax = 10, 
    k_per_logint = 20, 
    CAMBAccuracyBoost=1.1,
    CLAccuracyBoost = 1.0, 
    CLIntegrationAccuracy = 1):   
    
    init_cosmolike(external_nz_modeling,mod_nz)

    (log10k_interp_2D, z_interp_2D, lnPL, lnPNL, G_growth, z_interp_1D, chi) = get_camb_cosmology(omegam=omegam, 
        omegab=omegab, H0=H0, ns=ns, As_1e9=As_1e9, w=w, w0pwa=w0pwa, AccuracyBoost=AccuracyBoost, kmax=kmax,
        k_per_logint=k_per_logint, CAMBAccuracyBoost=CAMBAccuracyBoost)

    CLAccuracyBoost = CLAccuracyBoost * AccuracyBoost
    CLSamplingBoost = CLAccuracyBoost * AccuracyBoost
    CLIntegrationAccuracy = max(0, CLIntegrationAccuracy + 5*(AccuracyBoost-1.0))
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
    ci.set_nuisance_shear_calib(M = M)
    ci.set_nuisance_shear_photoz(bias = shear_photoz_bias)
    ci.set_nuisance_ia(A1 = A1, A2 = A2, B_TA = BTA)

    if baryon_sims is None:
        ci.reset_bary_struct()
    else:
        ci.init_baryons_contamination(sim = baryon_sims)
        
    (xip, xim) = ci.xi_pm_tomo()    
    return (ci.get_binning_real_space(), xip, xim)


# (theta0, xip0, xim0) = xi(external_nz_modeling=0)
# (theta1, xip1, xim1) = xi(external_nz_modeling=1)

# print(xip0[:,0,0]==xip1[:,0,0])

# plt.plot(theta0,theta0*xip0[:,0,0]*1e4)
# plt.plot(theta0,theta0*xip1[:,0,0]*1e4,ls="--")
# plt.xscale("log")
# plt.savefig("test.pdf")
