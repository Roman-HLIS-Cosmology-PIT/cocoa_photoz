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

# Libs to explore the derivative of xi_+ wrt n(z)
from scipy.signal import unit_impulse
from multiprocessing import Pool

path = "../../../external_modules/data/lsst_y1/lsst_y1_source.nz"
nz_fid = np.genfromtxt(path)
(theta_fid, xip_fid, xim_fid) = \
    cp.xi(external_nz_modeling=1,mod_nz=nz_fid)


epsilon = 0.01
tomo_bin = 1
n_z = nz_fid.shape[0]
z_vals = nz_fid[:,0]

def compute_derivative(args):
    z_idx = args
    nz_per = nz_fid.copy()
    nz_per[z_idx, tomo_bin] += epsilon
    nz_per[:,tomo_bin] /= np.trapz(y=nz_per[:,tomo_bin], x=z_vals)
    (theta_per, xip_per, xim_per) = \
        cp.xi(external_nz_modeling=1,mod_nz=nz_per)
    dxi_dn = (xip_per[:,0,0] - xip_fid[:,0,0]) / epsilon
    print(z_idx)
    return dxi_dn

jobs = list(range(n_z))

n_cores=4
with Pool(processes=n_cores) as pool:
    results = pool.map(compute_derivative, jobs)