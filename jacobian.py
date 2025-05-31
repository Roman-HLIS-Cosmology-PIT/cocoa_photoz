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
from math import factorial
# Libs to explore the derivative of xi_+ wrt n(z)
from scipy.signal import unit_impulse
from multiprocessing import Pool

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path_jacob = "/home/grads/data/Diogo/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz/test.txt" # where to save the jacobian matrix
path = "../../../external_modules/data/lsst_y1/lsst_y1_source.nz"
nz_fid = np.genfromtxt(path)
n_z = nz_fid.shape[0]
n_tomo = nz_fid.shape[1] - 1
z_vals = nz_fid[:,0]
combs = int(factorial(n_tomo+2-1)/(2*factorial(n_tomo-1))) # number of non-repeated combinations of i,j in n_tomo
jacob_dim1 = cp.ntheta * combs # size of the 2pt function
# jacob_dim2 = 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(theta_fid, xip_fid, xim_fid) = cp.xi(external_nz_modeling=1,mod_nz=nz_fid)

epsilon = 0.01

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_derivative(args):
    z_idx, tomo_z_bin = args
    nz_per = nz_fid.copy()
    nz_per[z_idx, tomo_z_bin+1] += epsilon
    nz_per[:,tomo_z_bin+1] /= np.trapz(y=nz_per[:,tomo_z_bin+1], x=z_vals)
    (theta_per, xip_per, xim_per) = cp.xi(external_nz_modeling=1,mod_nz=nz_per)
    
    dxi_dn = np.array([((xip_per[:,tbi,tbj] - xip_fid[:,tbi,tbj]) / epsilon) for tbi in range(n_tomo) for tbj in range(n_tomo) if tbj>=tbi]).reshape(1,jacob_dim1)

    return dxi_dn
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

jobs = [(zi, tb) for tb in range(n_tomo) for zi in range(n_z)]

derivs_map = map(compute_derivative,jobs)

with open(path_jacob,"a") as f:
    derivs = list(derivs_map)
    np.savetxt(f,derivs)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
