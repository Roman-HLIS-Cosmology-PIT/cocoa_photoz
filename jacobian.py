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
import cocoa_photoz as cp

from scipy.signal import unit_impulse
from multiprocessing import Pool

# path = "../../../external_modules/data/lsst_y1/lsst_y1_source.nz"
path = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/lsst_y1/data/lsst_y1_source.nz"
nz_fid = np.genfromtxt(path)
(theta_fid, xip_fid, xim_fid) = cp.xi(external_nz_modeling=1,mod_nz=nz_fid)

epsilon = 0.01
n_z = nz_fid.shape[0]
n_tomo = nz_fid.shape[1] - 1
z_vals = nz_fid[:,0]

def compute_derivative(args):
    z_idx, tomo_bin = args
    nz_per = nz_fid.copy()
    nz_per[z_idx, tomo_bin] += epsilon
    nz_per[:,tomo_bin] /= np.trapz(y=nz_per[:,tomo_bin], x=z_vals)
    (theta_per, xip_per, xim_per) = cp.xi(external_nz_modeling=1,mod_nz=nz_per)
    dxi_dn = (xip_per[:,0,0] - xip_fid[:,0,0]) / epsilon
    print("z_idx, tomo_bin: ",z_idx,tomo_bin)
    return dxi_dn

jobs = [(zi, tb) for tb in range(1,n_tomo+1) for zi in range(n_z)]

list(map(compute_derivative,jobs))



# with Pool() as pool:
#     results = pool.map(compute_derivative, jobs)