"""
Diogo H. F. de Souza - Fri Mar 28 2025
Get the mean of redshift distribution for 9 millions of Roman scenarios
"""

## LIBS ##
import numpy as np
import h5py

###################################
###### REDSHIFT DISTRIBUTION ######
###################################

nzs = h5py.File('roman_nz_realizations/sc1b_d4/SVSN/nz_samples__LHC0_pointZ_1e6_Roman_sc1b_d4.h5','r') # nz simulations: z and nzs
mean_nz_unorm = np.genfromtxt('results/zdistris/05_mean_roman_nz.txt') # unormalized

def ni_gal():
    """
    input: unormalized redshift distribution
    outputs: (used in CosmoCov inifile): source_n_gal = lens_n_gal = n_gal = fi/∑_{i=1}^9 (fi) where fi = ∫n_i(z) dz
    """
    int=[]
    for i in range(mean_nz_unorm.shape[1]):
        fi = [np.trapz(y=mean_nz_unorm[:,i], x=nzs["zbinsc"][:])]
        int+=fi

    tot = np.sum(int)
    fractions = int/tot
    return int,fractions # integral , ni_gal

mean_nz_norm = mean_nz_unorm.copy()
int , _ = ni_gal()
for i in range(mean_nz_norm.shape[1]):
    mean_nz_norm[:,i] = mean_nz_norm[:,i]/int[i]
np.savetxt('results/zdistris/05_mean_roman_nz_normalized_sc1bd4.txt',mean_nz_norm)

def make_z_nz_combined():
    zbinsc = nzs["zbinsc"][:] # z center (cocoa)
    zbinsmin = nzs["zbins"][:46] # z bin min (cosmocov) 

    zbinsc_reshaped = zbinsc.reshape(-1, 1)  # shape (46, 1)
    zbinsmin_reshaped = zbinsmin.reshape(-1, 1)  # shape (46, 1)

    combined_zbinsc = np.hstack((zbinsc_reshaped, mean_nz_norm))
    combined_zbinsmin = np.hstack((zbinsmin_reshaped, mean_nz_norm))

    np.savetxt("results/zdistris/mean_roman_nzc_normalized_sc1bd4.nz", combined_zbinsc) # z mean + nzs (cocoa)
    np.savetxt("results/zdistris/mean_roman_nzmin_normalized_sc1bd4.nz", combined_zbinsmin) # z min + nzs (cosmocov)
    return None
make_z_nz_combined()