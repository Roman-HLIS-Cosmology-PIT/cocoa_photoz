# Diogo Souza - Jul 14 2025
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import h5py
from sklearn.decomposition import PCA
from getdist import MCSamples, plots
from scipy.stats import norm
from numpy import linalg as la

from cobaya.yaml import yaml_load
from yaml import safe_load, dump
from cobaya.model import get_model
import sys, os
sys.path.insert(0, os.environ['ROOTDIR']+'/external_modules/code/CAMB/build/lib.linux-x86_64-'+os.environ['PYTHON_VERSION'])
import camb
from camb import model

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sim', type=str, required=True, help='Simulation name or tag')
args = parser.parse_args()

sim = args.sim
print(f"Simulation selected: {sim}")

plt.rcParams['figure.figsize'] = (3.5, 2.5)

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
matplotlib.rcParams['font.size'] = 15
##################################################################
##################################################################
##################################################################
colors1 = [
    "#1b5f6f",  # dark teal
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#882255",  # wine
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#999999",  # grey
    "#117733",  # dark green
]

camb_path = os.path.dirname(camb.__file__)
cocoa_path = '/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa' # Change for your path 
path = f'{cocoa_path}/cocoa_photoz'
print('Using CAMB %s installed at %s'%(camb.__version__,camb_path))
##################################################################
##################################################################
##################################################################
class SURVEY:
    def __init__(self,survey,use_r=None):
        self.survey = survey
        self.use_r = use_r

    # Global parameters
    def nz_params(self):
        if self.survey == 'des':
            Nz = 299 #len(z)
            Nt = 4
        elif self.survey == 'roman':   
            Nz = 46 #len(z)
            Nt = 9
        return Nz,Nt

    # Auxiliary functions
    def get_nzs(self):
        Nz,Nt = self.nz_params()
        if self.survey == 'des':
            nz = f'{path}/roman_nz_realizations/Fisher_matrix/Tz_realizations_WZ_bq_pile3_0d01.npy' # shared by Boyan Yin
            nz = np.load(nz)             ## shape = (10095, 4, 300) = (Ns, Nt, Nz) = (# of simulations, # of tomo bins, # of redshift)
            nz = nz[:,:,1:]              ## shape = (10095, 4, 299) | (first :) => all 10095 simulations, (second :) => all 4 tomo bins, (third ,1:) => all redshift bins except the first where nz is 0.0.
            Ns = np.shape(nz)[0]         ## Number of simulations: 10095
            nz = nz.reshape(Ns, Nt * Nz) ## shape = (10095, 1196) = (Ns, Nt*Nz)
            min_z   = 0.01
            max_z   = 3
            delta_z = (max_z - min_z) / Nz
            zbins   = np.arange(min_z,max_z+delta_z,delta_z)
            z  = zbins[:-1]+(zbins[1]-zbins[0])/2.
            return z, nz

        elif self.survey == 'roman':
            sc = 'sc1bd4'
            nz = h5py.File(f'{path}/roman_nz_realizations/{sc}/nz_samples_LHC0_pointZ_1e6_Roman_{sc}.h5','r') # shared by Boyan Yin
            z = np.array(nz['zbinsc'])
            return z, nz

    def roman_nbar(self):    
        nbar = f'{path}/n_nbar_ndiff_covndiff/nbar_sc1b_d4.nz'
        return nbar

    def cov_of_ndiff(self):
        Nz,Nt = self.nz_params()
        z,nz = self.get_nzs()

        if self.survey == 'roman':
            # Never use r=10^6 on login node.
            if self.use_r is not None:
                #np.random.seed(42)
                #sims_indexes = np.random.randint(1, 10**6 + 1, size=self.use_r)
                # tot = len(sims_indexes)
                sims_indexes = self.use_r

                rows = []
                print(f"Total # of realizations r = {self.use_r}\nNormalizing nz for r:")
                for jj, j in enumerate(sims_indexes):
                    row = []
                    for i in range(Nt):
                        norm = np.trapz(y=nz[f'bin{i}'][j], x=z)
                        normalized = nz[f'bin{i}'][j] / norm
                        row.append(normalized)
                    #row shape after appending all bins: (Nt,Nz)
                    print(jj, end=' ')
                    rows.append(np.hstack(row))
                #rows shape after appending all realizations: (tot,Nt*Nz)
                nz = np.vstack(rows) # (tot,Nt*Nz)
                nbar = np.mean(nz,axis=0)
                ndiff = nz - nbar
                Cn = np.einsum('ij,ik->jk',ndiff,ndiff) / (ndiff.shape[0]-1) # equivalent to Cn = ndiff.T @ ndiff / sample size - 1
            else:
                # Use Cn from all 10^6
                # Cn = np.load(f'{path}/n_nbar_ndiff_covndiff/cov_ndiff_sc1b_d4.npy') #TODO: path changed    
                print("USING COVARIANCE OF REALIZATIONS:")
                path_to_Cn = f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/sc1bd4/Cn_sc1bd4.txt'
                print(path_to_Cn)
                Cn = np.genfromtxt(path_to_Cn) # TODO: path changed - !!!!  AUTOMATIZAR !!!!   TODO: INTEGRAR COM OUTPUTS DE PIP_1_NNORM_NBAR_NDIFF_CN.PY
                ndiff = None # Because we can't load (414,1M) data here
        elif self.survey == 'des':
            ### nzs for DES already are normalized
            ### set parameter below to True if you want to check
            check_des_nz_normalization = False
            if check_des_nz_normalization:
                rows = []
                for j in range(10095):
                    row = []
                    for i in range(Nt):
                        start, end = Nz*i, Nz*(i+1)
                        norm = np.trapz(y=nz[j,start:end], x=z)
                        normalized = nz[j,start:end] / norm
                        row.append(normalized)
                    rows.append(np.hstack(row))
                nz = np.vstack(rows)

            nbar = np.mean(nz,axis=0)
            ndiff = nz - nbar
            Cn = np.einsum('ij,ik->jk',ndiff,ndiff) / (ndiff.shape[0]-1) # equivalent to Cn = ndiff.T @ ndiff / sample size - 1
        self.Cn = Cn
        self.ndiff = ndiff    
        return self.Cn, self.ndiff
    
    # Unweighted PCA    
    def pcs(self,method='eig'):
        Nz,Nt = self.nz_params()
        Nd = Nz*Nt
        if method == 'eig':
            eig_vals, eig_vecs = np.linalg.eig(self.Cn)
            eigvals, eigvecs = np.linalg.eig(self.Cn)
        elif method == 'eigh':
            eigvals, eigvecs = np.linalg.eigh(self.Cn)
        elif method == 'svd':    
            U, S, VT = np.linalg.svd(self.Cn)
            eigvecs, eigvals = U,S
        elif method == 'sklearn':
            pca = PCA()
            pca.fit(self.Cn)
            eigvecs, eigvals = pca.components_, pca.singular_values_
            explained_variance = pca.explained_variance_
            explained_variance_ratio = pca.explained_variance_ratio_
        ######
        if method != 'sklearn':
            eigvals = np.real(eigvals)
            eigvecs = np.real(eigvecs)
            idx     = np.argsort(eigvals)[::-1] # 1. Get indices to sort eigenvalues in descending order
            eigvals = eigvals[idx]              # 2. Sort eigenvalues (largest to smallest)
            eigvals = np.maximum(0,eigvals)     # 3. Clamp any small negative eigenvalues to 0 (often due to numerical errors)
            eigvecs = eigvecs[:, idx]           # 4. Reorder columns of eigenvectors to match sorted eigenvalues
            return eigvecs, eigvals
        elif method == 'sklearn':
            return eigvecs.T, eigvals, explained_variance,explained_variance_ratio
    
    def forward_pca(self, params_values,nbar_path,pcs_path,npcs_nz):
        nbar = np.genfromtxt(nbar_path)
        U = np.genfromtxt(pcs_path)[:,:npcs_nz]
        s = nbar[:,1:].shape
        z = nbar[:,0]
        # Model: n(z) = <n>(z) + α_1*PC_1(z) + α_2*PC_2(z) + ... + α_n*PC_n(z)
        if npcs_nz > 0:
            alphas = np.array([params_values.get("roman_alpha_"+str(i+1)) for i in range(npcs_nz)])
            correction = (alphas * U).sum(axis=1)
    
            nz_model = nbar[:,1:].T.flatten() + correction
            nz_model = nz_model.reshape(s[::-1]).T
            nz_model = np.column_stack((z,nz_model))
            return nz_model
        else:
            return nbar

    # Credit: Gary M. Bernstein
    def nearestPD(self,A):
        """Find the nearest positive-definite matrix to input
        A Python/Numpy port of John D'Errico's `nearestSPD` MATLAB code [1], which
        credits [2].
        [1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd
        [2] N.J. Higham, "Computing a nearest symmetric positive semidefinite
        matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6
        """
        B = (A + A.T) / 2
        _, s, V = np.linalg.svd(B)
        H = np.dot(V.T, np.dot(np.diag(s), V))
        A2 = (B + H) / 2
        A3 = (A2 + A2.T) / 2
        if self.isPD(A3):
            return A3
        spacing = np.spacing(np.linalg.norm(A))
        I = np.eye(A.shape[0])
        k = 1
        while not self.isPD(A3):
            mineig = np.min(np.real(np.linalg.eigvals(A3)))
            A3 += I * (-mineig * k**2 + spacing)
            k += 1
        return A3

    def isPD(self,B):
        """Returns true when input is positive-definite, via Cholesky"""
        try:
            _ = np.linalg.cholesky(B)
            return True
        except np.linalg.LinAlgError:
            return False              
    
    def getModes(self,D, Cn, chisq_threshold=0):
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

        You can safely use fewer than `nModes` of the returned X and U if you 
        think the discarded `dchisq` values are small enough.'''

        thr = 1e-25     # Dynamic range threshold for eigenvalues

        # Symmetrize and sqrt D
        D = 0.5*(D + D.T)
        Dval, Dvec = np.linalg.eigh(D)
        # Any negatives are numerical problems
        Dval = np.maximum(0., Dval)

        # And Cn
        Cnval, Cnvec = np.linalg.eigh(Cn)
        Cnval = np.maximum(0., Cnval)

        # Build the keystone matrix and its SVD
        M = np.einsum('i,ji,jk,k->ik',np.sqrt(Cnval), Cnvec, Dvec, np.sqrt(Dval))
        Um,s,Vmt = np.linalg.svd(M)

        # Sort SV's and throw away unwanted ones
        order = np.argsort(s*s)  # increasing order
        kill = np.count_nonzero(np.cumsum(s[order]**2) < chisq_threshold)
                                            
        resid = np.sum((s[order[:kill]]**2))

        keep = np.flip(order[kill:])  # Decreasing influence order
        Um = Um[:,keep]
        Vm = Vmt.T[:,keep]
        s = s[keep]

        # Build the encoder
        tmp = np.where(Cnval>thr*np.max(Cnval), 1/np.sqrt(Cnval), 0.)
        X = np.einsum('ij,i,ki->jk',Um, tmp, Cnvec)  # X is (M,N)
        # And the decoder
        tmp = np.where(Dval>thr*np.max(Dval), 1/np.sqrt(Dval), 0.)
        U = np.einsum('ji,i,ik,k->jk',Dvec, tmp, Vm,s) # Vm s
        
        return X, U, s*s, resid, Dvec, Dval, Cnvec, Cnval, Um, Vm, M
    
    ###############################
    ############ CoCoA ############
    ###############################
    
##################################################################
##################################################################
##################################################################
roman_obj = SURVEY('roman',use_r=None)
Nz,Nt     = roman_obj.nz_params()
z,nz      = roman_obj.get_nzs()
nbar_file = roman_obj.roman_nbar()
Cn, ndiff = roman_obj.cov_of_ndiff() # must call before pcs

eig_vecs, eig_vals = roman_obj.pcs('eig')
eih_vecs, eih_vals = roman_obj.pcs('eigh')
svd_vecs, svd_vals = roman_obj.pcs('svd')
# np.savetxt('U.txt',eih_vecs) # Jul 18 using on MCMC

skl_vecs, skl_vals,\
explained_variance, explained_variance_ratio = roman_obj.pcs('sklearn')

def fisher(weights=True):
    if weights:
        stencil_5pt   = f'{cocoa_path}/jacobian_outs/stencil_5pt'
        step_sizes    = '1.000000e-03'
        relative_path = 'stencil_5pt'
        inv_cov       = np.genfromtxt(f'{stencil_5pt}/inv_cov.txt')
        deriv_ξpm     = np.genfromtxt(f'{cocoa_path}/jacobian_outs/{relative_path}/deriv_ξpm_{step_sizes}.txt')
        fisher_ξpm    = deriv_ξpm@inv_cov[:deriv_ξpm.shape[1],:deriv_ξpm.shape[1]]@deriv_ξpm.T
        print("USING WEIGHTS")
    if not weights:
        print("NOT USING WEIGHTS")
        fisher_ξpm = np.identity(Cn.shape[0])   

    if not roman_obj.isPD(fisher_ξpm):
        print('fixing cov')
        F=roman_obj.nearestPD(fisher_ξpm)
    else:
        F=fisher_ξpm 
    return F
##################################################################
##################################################################
##################################################################
npcs_tot  = 414
npcs_used = 414

print("--------------------------------------------------")
print("-------------------nz SIGNATURE-------------------")
print(f"realization # {sim}")
print("-------------------nz SIGNATURE-------------------")
print("--------------------------------------------------")

nbar_path     = f'{cocoa_path}/projects/roman_real/data/sc1bd4/nbar_sc1bd4.nz'
nz_fid_path   = f'{cocoa_path}/projects/roman_real/data/sc1bd4/nz_sc1bd4r{sim}.nz'

nz_fid     = np.genfromtxt(nz_fid_path)
nbar       = np.genfromtxt(nbar_path)
z          = nz_fid[:,0]
ndiff      = nz_fid[:,1:] - nbar[:,1:] # (46,9)
##################################################################
##################################################################
##################################################################
def projection(encoding,decoding):
    decoding=decoding[:,:npcs_used]
    # params_values={}
    # for i in range(npcs_used):
    #     alpi = np.dot(ndiff.T.flatten(),decoding[:,i])
    #     # print(f'  roman_alpha_{i+1}: {alpi}')
    #     params_values.update({f"roman_alpha_{i+1}": alpi})

    params_values={}
    u = ndiff.T.flatten() @ encoding.T # Eq 21 of 2506.00758 (compression)
    for i in range(npcs_used):
        print(f'  roman_alpha_{i+1}: {u[i]}')
        params_values.update({f"roman_alpha_{i+1}": u[i]})

    alphas = np.array([params_values.get("roman_alpha_"+str(i+1)) for i in range(npcs_used)])
    correction = (alphas * decoding).sum(axis=1)

    s = nbar[:,1:].shape
    z = nbar[:,0]
    nz_model = nbar[:,1:].T.flatten() + correction # Eq 22 of 2506.00758 (decompression)
    nz_model = nz_model.reshape(s[::-1]).T
    nz_model = np.column_stack((z,nz_model))

    nz_model[nz_model < 0] = 0
    norms = np.trapz(y=nz_model[:,1:],x=z,axis=0)
    nz_model[:,1:] /= norms    

    for i in range(1,10):
        plt.plot(z,nz_fid[:,i],c="C0",label="fid" if i==1 else None)
        plt.plot(z,nz_model[:,i],c="k",ls="--",label="recon" if i==1 else None)
        plt.plot(z,nbar[:,i],c="orange",ls=":",label="mean" if i==1 else None)

    plt.legend(loc="upper right")
    plt.savefig(f"test{npcs_used}.pdf")

    return None

def execute_end_to_end_pipeline(extra_info=False):
    F = fisher(weights=True)
    E, D, ss, resid, Dvec, Dval, Cnvec, Cnval, Um, Vm, M = roman_obj.getModes(F, Cn, chisq_threshold=0)
    projection(encoding=E,decoding=D)

    if extra_info:
        idx_desc = np.argsort(Cnval)[::-1]
        Cnval = Cnval[idx_desc]
        Cnvec = Cnvec[:, idx_desc]

        X = D@E.T
        ss_tot = np.sum(ss)
        test = ss_tot - np.cumsum(ss)
        np.savetxt(f"E_{sim}.txt",E)
        np.savetxt(f"D_{sim}.txt",D)
        np.savetxt(f"X_{sim}.txt",X)
        np.savetxt(f"Cnvec_{sim}.txt",Cnvec)
        np.savetxt(f"test_{sim}.txt",test)
    return None

# execute_end_to_end_pipeline(extra_info=True)

##################################################################
##################################################################
##################################################################

def chi2_vs_npcs():
    # Load the yaml as Python dict
    # del info_txt, info, point
    input_file  = f"r{sim}_weights0.yaml"
    output_file = f"chi2_vs_npcs_r{sim}_weights0.txt"
    info_txt    = f'{path}/yamls/roman_pca/chi2_vs_npcs_weighted/{input_file}'
    print("----------------------------------")
    print(f"[INPUT] USING YAML {info_txt}")
    print(f"[OUTPUT] SAVING AT {output_file}")
    print("----------------------------------\n")
    
    with open(info_txt) as f:
        info = yaml_load(f)

    # Adjust the paths for theory (camb) and likelihood (roman) codes
    del info["theory"]["camb"]["path"]
    info["packages_path"] = camb_path
    info["likelihood"]["roman_real.roman_real_cosmic_shear"]["path"] = f'{cocoa_path}/external_modules/data/roman_real/sc1bd4'

    open(f"{cocoa_path}/chi2_vs_npcs_weighted/{output_file}","w").close()

    with open(f"{cocoa_path}/chi2_vs_npcs_weighted/{output_file}","a") as f:
        f.write("# PCs, χ²\n")
        f.write(f"# Using yaml {info_txt}\n")
        for i in range(npcs_used):
            info["likelihood"]["roman_real.roman_real_cosmic_shear"]["npcs_nz"] = i
            model = get_model(info)
            varied_points = {} # The yaml "example_unweighted.yaml" do not contain any varied parameter! 
            loglikes_i = model.loglikes(varied_points,as_dict=False,return_derived=False)[0]
            chi2_i = -2*loglikes_i
            print(f"# of PCs: {info['likelihood']['roman_real.roman_real_cosmic_shear']['npcs_nz']}, χ²: {chi2_i}\n",end='')
            f.write(f"{i}, {chi2_i}\n")
    return None        

# chi2_vs_npcs()