"""
Diogo H. F. de Souza - Thu Feb 27 2025
Modification to the original code `nz_compression.py`
Original code credit: Troxel, Boyan, Clare, ...
"""

## LIBS ##
# import fitsio as fio
import numpy as np
from numpy import linalg as la
from numpy.lib.recfunctions import stack_arrays
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import pylab
from scipy.special import softmax
import h5py
# import twopoint
import os

# SURVEY = 'ROMAN' # or DES
SURVEY = 'DES' # or ROMAN

################################
###### Fisher matrix ######
################################

def fisher(survey):
    if survey == 'DES': 
        # fisher = 'results/fisher_matrix/des_y3/fisher.txt'
        fisher = 'fisher2_all_with_scale_cuts.txt'

    elif survey == 'ROMAN':
        fisher = 'results/fisher_matrix/roman_real/fisher.txt'

    D = np.genfromtxt(fisher)
    return D    

def get_nzs(survey):
    if survey == 'DES':
        ## nz simulations shared by Boyan
        nz_file_DES = 'roman_nz_realizations/Fisher_matrix/Tz_realizations_WZ_bq_pile3_0d01.npy'
            
        nzs = np.load(nz_file_DES) ## shape = (10095, 4, 300) = (Ns, Nt, Nz)
        nzs = nzs[:,:,1:]          ## shape = (10095, 4, 299) DHFS: For all 10095 simulations (first :) and for all 4 tomo bins (second :) eliminate eliminate the first element (third ,1:) because the first element of all nzs is 0.0.
        Ns = np.shape(nzs)[0]      ## Number of simulations: 10095
        Nt = np.shape(nzs)[1]      ## Number of tomographic bins: 4
        Nz = np.shape(nzs)[2]      ## Number of redshifts: 299
        nzs_DES = nzs.reshape(Ns, Nt * Nz)
        return nzs_DES

    elif survey == 'ROMAN':
        # open simulation: this file contains z and nzs
        nz_file_ROMAN = 'roman_nz_realizations/sc1b_d4/SVSN/nz_samples__LHC0_pointZ_1e6_Roman_sc1b_d4.h5'
        nzs_ROMAN_ = h5py.File(nz_file_ROMAN,'r') ## nzs in dict format
        sim_j = lambda j: [nzs_ROMAN_[f'bin{i}'][j,:] for i in np.arange(9)] ## The jth simulation for all 9 tomographic bins.
        # Ns_Roman = 10095 ## DHFS - Taken Ns for Roman to be equal to Ns for DES - disclaimer: should use all 9M simulations for Roman
        Ns_Roman = 50000 ## DHFS - Taken Ns for Roman to be equal to Ns for DES - disclaimer: should use all 9M simulations for Roman
        nzs_ROMAN = np.stack(([sim_j(j) for j in np.arange(Ns_Roman)])) ## nzs in array format
        nzs_ROMAN = nzs_ROMAN.reshape(Ns_Roman, 9 * 46)
        return nzs_ROMAN

########################################################
##### TEST 
########################################################
# nz_file_DES = '/gpfs/projects/MirandaGroup/Diogo/ROMAN-NZ-PROJECT/photoz_uz/Tz_realizations_WZ_bq_pile3_0d01.npy'  ## DHFS: Shared by Boyan, Troxel - this file is for DES
# nz_file_ROMAN = '/gpfs/projects/MirandaGroup/Diogo/ROMAN-NZ-PROJECT/sc1b_d4/SVSN/nz_samples__LHC0_pointZ_1e6_Roman_sc1b_d4.h5' # open simulation: this file contains z and nzs

# nzs_DES   = np.load(nz_file_DES)
# nzs_ROMAN_ = h5py.File(nz_file_ROMAN,'r') ## nzs in dict format

# sim_j = lambda j: [nzs_ROMAN_[f'bin{i}'][j,:] for i in np.arange(9)] ## The jth simulation for all 9 tomographic bins.
# Ns_Roman = 2 ## DHFS - Taken Ns for Roman to be equal to Ns for DES - disclaimer: should use all 9M simulations for Roman
# nzs_ROMAN = np.stack(([sim_j(j) for j in np.arange(Ns_Roman)])) ## nzs in array format

# for j in np.arange(Ns_Roman):
#     for i in np.arange(9):
#         plt.plot(nzs_ROMAN_['zbinsc'][:],nzs_ROMAN[j,i,:])
# plt.savefig('./test.pdf')


# N = (3-0)/0.01
# z = np.linspace(0,3,int(N)-1)

# for j in np.arange(5):
#     for i in np.arange(4):
#         plt.plot(z,nzs[j,i,:])
# plt.savefig('./test.pdf')
# print(x1==x2)
# print(x1)
# print(x2)
# print(nzs.shape)
# nzs = nzs.reshape(np.shape(nzs)[0], 4*299)
# print(nzs.shape)

########################################################
##### TEST 
########################################################

def nearestPD(A):
    """Find the nearest positive-definite matrix to input
    A Python/Numpy port of John D'Errico's `nearestSPD` MATLAB code [1], which
    credits [2].
    [1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd
    [2] N.J. Higham, "Computing a nearest symmetric positive semidefinite
    matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6
    """
    B = (A + A.T) / 2
    _, s, V = la.svd(B)
    H = np.dot(V.T, np.dot(np.diag(s), V))
    A2 = (B + H) / 2
    A3 = (A2 + A2.T) / 2
    if isPD(A3):
        return A3
    spacing = np.spacing(la.norm(A))
    I = np.eye(A.shape[0])
    k = 1
    while not isPD(A3):
        mineig = np.min(np.real(la.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1
    return A3


def isPD(B):
    """Returns true when input is positive-definite, via Cholesky"""
    try:
        _ = la.cholesky(B)
        return True
    except la.LinAlgError:
        return False


def getModes(D, Cn, chisq_threshold=0.1):
    '''Calculate the compression matrix from n to u and the
    modes U that multiply each U.  
    Parameters:
    `D`:      The matrix J^T C_c^{-1} J, where C_c is the covariance
              matrix of the observables c, and J is the Jacobian dc/dn,
              shape=(N,N)
    `Cn`:     Covariance matrix of the n(z) parameter vector n, shape=(N,N)
    `chisq_threshold`:  Your desired upper bound for the mean error in chisq
              resulting from the compression.
    Returns:
    `X`:      The compression matrix, u = X @ n, shape=(nModes,N).  The
              covariance matrix of the u's will be the identity.
    `U`:      The decoding matrix, n' = U @ u, shape=(N,nModes)
    `dchisq`: The mean chisq shift generated by each mode, shape=(nModes)
    `resid`:  Mean amount of chisq shift caused by compression.
    ######### DHFS START #########
    Returns:
    `X=E`:    The Encoding matrix, u = E @ n, shape=(nModes,N).  The
              covariance matrix of the u's will be the identity.
    `U=D`:    The Decoding matrix, n' = D @ u, shape=(N,`nModes`)
    ######### DHFS END #########

    You can safely use fewer than `nModes` of the returned X and U if you 
    think the discarded `dchisq` values are small enough.'''

    thr = 1e-25     # Dynamic range threshold for eigenvalues

    ## DHFS: Fisher matrix J^T C_c^{-1} J where J = dc/dn
    # Symmetrize and sqrt D
    D = 0.5*(D + D.T)
    Dval, Dvec = np.linalg.eigh(D)
    # Any negatives are numerical problems
    Dval = np.maximum(0., Dval)

    ## DHFS: Covariance matrix of n(z)
    # And Cn
    Cnval, Cnvec = np.linalg.eigh(Cn)
    Cnval = np.maximum(0., Cnval)

    # Build the keystone matrix and its SVD
    M = np.einsum('i,ji,jk,k->ik',np.sqrt(Cnval), Cnvec, Dvec, np.sqrt(Dval))
    Um,s,Vmt = np.linalg.svd(M)

    # Sort SV's and throw away unwanted ones
    order = np.argsort(s*s)  # increasing order
    kill = np.count_nonzero(np.cumsum(s[order]**2) < chisq_threshold)
                                           
    resid = np.sum((s[order[:kill]]**2)) ## DHFS: Gary's Eq. 15 - have to confirm!

    keep = np.flip(order[kill:])  # Decreasing influence order
    Um = Um[:,keep]
    Vm = Vmt.T[:,keep]
    s = s[keep]

    # Build the encoder
    tmp = np.where(Cnval>thr*np.max(Cnval), 1/np.sqrt(Cnval), 0.)
    X = np.einsum('ij,i,ki->jk',Um, tmp, Cnvec)  # X is (M,N)
    # And the decoder
    tmp = np.where(Dval>thr*np.max(Dval), 1/np.sqrt(Dval), 0.)
    U = np.einsum('ji,i,ik,k->jk',Dvec, tmp, Vm,s)
    
    return X, U, s*s, resid

# Load nzs and calculate mean, deviation from mean, and covariance
chisq_threshold=0.15
# chisq_threshold=1e-25

n=get_nzs( SURVEY )
nbar = np.mean(n, axis=0)
ndiff = n - nbar
Cn = np.einsum('ij,ik->jk',ndiff,ndiff) / n.shape[0]
print('Cn shape:',Cn.shape)

D = fisher( SURVEY )

# Load the D matrix (fisher matrix of n(z_i) variations) - its eigenvalues drop rapidly too
#D = np.genfromtxt('chain_nzsample_source_unsmooth_coarse_max3.txt')
# This has never turned out to be positive definite, needs to be 'corrected' first

if not isPD(D):
    print('fixing cov')
    D=nearestPD(D)

# Eigen* of fisher matrix D
print('D shape:',D.shape)

X,U,dchisq,resids = getModes(D, Cn, chisq_threshold=chisq_threshold)
print('Chisq kept:',dchisq,'discarded:',resids)
# print(U)

## Save eigenvectors/basis/modes to file
print('U shape 1:',U.shape)
if SURVEY=="DES":
    U = np.reshape(U.T, (np.shape(U.T)[0], 4, -1))
elif SURVEY=="ROMAN":
    U = np.reshape(U.T, (np.shape(U.T)[0], 9, -1))
# print('U:',U)
print('U shape 2:',U.shape)
np.savez('./U_source.npz', U=U, perbin=0)

# Encode
u = ndiff @ X.T   # u is (Nz,M)
nEig = np.shape(u)[1]
print('nEig: ', nEig)

chisq_kept = np.array([dchisq])
np.savetxt('chisq_kept.txt', chisq_kept)
chisq_discard = np.array([resids])
np.savetxt('chisq_discard.txt', chisq_discard)