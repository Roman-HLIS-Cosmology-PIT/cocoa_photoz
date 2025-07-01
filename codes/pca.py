import numpy as np
import matplotlib.pyplot as plt
import h5py
from sklearn.decomposition import PCA
from getdist import MCSamples, plots

# Colors
colors = [
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

path='/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/'

nz = f'{path}/roman_nz_realizations/sc1b_d4/nz_samples_LHC0_pointZ_1e6_Roman_sc1b_d4.h5'
nz = h5py.File(nz,'r') 
z = np.array(nz['zbinsc'])
Nz = 46 #len(z)
Nt = 9
Nd = 414

def pcs(method='eig'):
    Cn = np.load(f'{path}/Cn_roman_sc1bd4.npy')
    Cn = 0.5*(Cn + Cn.T)
    if method == 'eig':
        eigvals, eigvecs = np.linalg.eig(Cn)
    elif method == 'eigh':
        eigvals, eigvecs = np.linalg.eigh(Cn)
    elif method == 'svd':    
        U, S, VT = np.linalg.svd(Cn)
        eigvecs, eigvals = U,S
    elif method == 'sklearn':
        pca = PCA()
        pca.fit(Cn)
        eigvecs, eigvals = pca.components_, pca.singular_values_
        explained_variance = pca.explained_variance_
        explained_variance_ratio = pca.explained_variance_ratio_
    ######
    if method != 'sklearn':
        eigvecs = np.real(eigvecs)
        eigvals = np.real(eigvals)
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]
        eigvecs = eigvecs[:, : Nd]
        eigvals = eigvals[: Nd]
        eigvals = np.maximum(0,eigvals)
        return eigvecs, eigvals
    elif method == 'sklearn':
        return eigvecs, eigvals, explained_variance,explained_variance_ratio 

eig_vecs, eig_vals = pcs('eig')
eigh_vecs, eigh_vals = pcs('eigh')
svd_vecs, svd_vals = pcs('svd')

eig_vals = eig_vals**2
eigh_vals = eigh_vals**2
svd_vals = svd_vals**2

skl_vecs, skl_vals,\
explained_variance, explained_variance_ratio = pcs('sklearn')

def plot_eigvalues_and_related_quantities(plot_type):
    plt.figure()
    if plot_type=='eig_vals':
        plt.plot(list(range(Nd)),eig_vals,'C0')
        plt.plot(list(range(Nd)),eigh_vals,'k--')
        plt.plot(list(range(Nd)),svd_vals,'lightgreen',ls=':')
        plt.plot(list(range(Nd)),skl_vals,'r-.')
        plt.xscale('log')
    elif plot_type == 'explained_variance_ratio':
        Ms = list(range(Nd))
        evr_eig = np.cumsum(eig_vals/np.sum(eig_vals))
        evr_eigh = np.cumsum(eigh_vals/np.sum(eigh_vals))
        evr_svd = np.cumsum(svd_vals/np.sum(svd_vals))
        evr_skl = np.cumsum(explained_variance_ratio)
        plt.plot(Ms,evr_eig,c=colors[0],ls='-',label='linalg.eig',lw=5)
        plt.plot(Ms,evr_eigh,c=colors[1],ls='--',label='linalg.eigh',lw=5)
        plt.plot(Ms,evr_svd,c=colors[2],ls='-.',label='linalg.svd',lw=5)
        plt.plot(Ms,evr_skl,c=colors[3],ls=':',label='sklearn',lw=5)
        plt.legend(loc='best',fontsize=15)
        plt.xlabel(r'$\text{Number of PCs:}~M$',fontsize=15)
        plt.ylabel(r'$\text{Explained variance ratio:}~\sum_{i=0}^{M}[\sigma_i^2/\sum_{i=0}^{414}\sigma_i^2]$',fontsize=15)
        plt.xticks(range(0,Nd + 1, 100))
        plt.tight_layout()
    elif plot_type == 'explained_variance_ratio_paper':
        M = 23 # maxium z
        Ms = list(range(Nd))[:M]
        x = np.cumsum(eig_vals/np.sum(eig_vals))[:M]*100
        plt.plot(Ms,x,c=colors[0],marker='o')
        plt.xlabel(r'$\text{Number of PCs:}~M$',fontsize=15)
        plt.ylabel(r'$\text{Explained variance ratio}~(\%)$',fontsize=15)
        plt.xlim(0,M)
        plt.ylim(np.min(x),np.max(x))
        plt.xticks(range(0,M + 1, 2))
        plt.tight_layout()
    plt.savefig('test.pdf')

# plot_eigvalues_and_related_quantities('explained_variance_ratio')

def plot_eigenvectors_singular_vectors(plot_type):
    k=4
    plt.figure()
    if plot_type=='test':
        plt.plot(z,eig_vecs[Nz*k:Nz*(k+1):,0],c=colors[0],ls='-',label='linalg.eig',lw=5)
        plt.plot(z,eigh_vecs[Nz*k:Nz*(k+1):,0],c=colors[1],ls='--',label='linalg.eigh',lw=5)
        plt.plot(z,svd_vecs[Nz*k:Nz*(k+1):,0],c=colors[2],ls='-.',label='linalg.svd',lw=5)
        plt.plot(z,skl_vecs.T[Nz*k:Nz*(k+1):,0],c=colors[3],ls=':',label='sklearn',lw=5)
        plt.legend(loc='best',fontsize=15)
        plt.xlabel(r'$z$',fontsize=15)
        plt.ylabel(r'$\mathrm{PC}_1^{\text{tomo}=1}$',fontsize=15)
        plt.xlim(np.min(z),np.max(z))
        plt.tight_layout()
    if plot_type == 'paper':
        plt.plot(z,eigh_vecs[Nz*k:Nz*(k+1):,0],c=colors[0],ls='-',lw=5,label=r'$\text{PC}_1$')
        plt.plot(z,eigh_vecs[Nz*k:Nz*(k+1):,1],c=colors[1],ls='--',lw=5,label=r'$\text{PC}_2$')
        plt.plot(z,eigh_vecs[Nz*k:Nz*(k+1):,2],c=colors[2],ls=':',lw=5,label=r'$\text{PC}_3$')
        # plt.plot(z,eigh_vecs[Nz*k:Nz*(k+1):,3],c=colors[3],ls='-.',lw=5)
        plt.xlabel(r'$z$',fontsize=15)
        plt.ylabel(r'$\mathrm{PC}s^{\text{tomo}=4}$',fontsize=15)
        plt.xlim(np.min(z),np.max(z))
        plt.legend(loc='best',fontsize=15)
        plt.tight_layout()
    plt.savefig('test.pdf')
    return None

# plot_eigenvectors_singular_vectors('paper')

def weights():
    alphas = [] # little us in the paper
    U, _ = pcs('eig') # (414=9*46, 414=#pcs)
    ndiff = np.load(f'{path}/ndiff_roman_sc1bd4.npy') # (1M,414=9*46)
    for t in range(9):
        n_sample = ndiff[0,Nz*t:Nz*(t+1)]
        PC = U[Nz*t:Nz*(t+1),0]
        alpha = np.dot(n_sample,PC.T) # projection <n,PC>
        alphas.append(alpha)
    return alphas


def plot_alpha_samples():
    M = 414
    Nsample = 1000

    U,_ = pcs('eig') # (414=9*46, 414=#pcs)
    ndiff = np.load(f'{path}/ndiff_roman_sc1bd4.npy') # shape (1M, 414=9*46)
    alphas = ndiff[:Nsample,:] @ U # Projection onto PCs: <n,PC>

    print(alphas.shape)

    names = ["alpha%s" %i for i in range(M)]
    labels = [fr"\alpha_{{{i+1}}}" for i in range(M)]

    chains = MCSamples(samples=alphas,names=names,labels=labels)

    g = plots.get_subplot_plotter()
    g.plots_1d(chains,["alpha%s" %i for i in range(5)],share_y=True)
    g.export('test.pdf')
    return None

# plot_alpha_samples()