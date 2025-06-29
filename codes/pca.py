import numpy as np
from numpy.linalg import eig
import matplotlib.pyplot as plt
import h5py
from getdist import MCSamples, plots

path='/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/'

nz = f'{path}/roman_nz_realizations/sc1b_d4/SVSN/nz_samples__LHC0_pointZ_1e6_Roman_sc1b_d4.h5'
nz = h5py.File(nz,'r') 
z = np.array(nz['zbinsc'])
Nz = 46 #len(z)
Nt = 9
Nd = 414

def pcs():
    d_cov = np.load(f'{path}/Cn_roman_sc1bd4.npy')
    eigvals, eigvecs = eig(d_cov)
    eigvecs = np.real(eigvecs)
    eigvals = np.real(eigvals)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    eigvecs = eigvecs[:, : Nd]
    eigvals = eigvals[: Nd]
    return eigvecs, eigvals

def weights():
    alphas = []
    U, _ = pcs() # (414=9*46, 414=#pcs)
    ndiff = np.load(f'{path}/ndiff_roman_sc1bd4.npy') # (1M,414=9*46)
    for t in range(9):
        n_sample = ndiff[0,Nz*t:Nz*(t+1)]
        PC = U[Nz*t:Nz*(t+1),0]
        alpha = np.dot(n_sample,PC.T)
        alphas.append(alpha)
    return alphas

alphas_1 = []
alphas_2 = []
# alphas_i = []
tomo = 0

M = 414
Nsample = 1000
# alphas_test = np.zeros((M, Nsample))  # rows: PCA modes, columns: Nsample

def alpha_sample():
    U,_ = pcs() # (414=9*46, 414=#pcs)
    ndiff = np.load(f'{path}/ndiff_roman_sc1bd4.npy') # shape (1M, 414=9*46)
    alphas = ndiff[:Nsample,:] @ U
    return alphas

alphas =  alpha_sample()

print(alphas.shape)
names = ["alpha%s" %i for i in range(M)]
labels = [fr"\alpha_{{{i+1}}}" for i in range(M)]

chains = MCSamples(samples=alphas,names=names,labels=labels)

g = plots.get_subplot_plotter()
g.plots_1d(chains,["alpha%s" %i for i in range(5)],share_y=True)
g.export('test.pdf')

###

# for i in range(9):
#     plt.figure()
#     plt.plot(z,U[46*i:46*(i+1),0],c='C0')
#     plt.plot(z,U[46*i:46*(i+1),1],c='k')
#     plt.plot(z,U[46*i:46*(i+1),2],c='r')
#     plt.savefig(f'test{i}.pdf')

# x = range(len(s))
# thr = [0.001,0.005,0.01]
# chisq_kept2 = np.genfromtxt(f'{path}/chisq_kept_{thr[1]}.txt')
# print(chisq_kept2)
# print(s[:6])
# for i in range(len(chisq_kept2)):
#     print(s[:6][i]/chisq_kept2[i])
# plt.plot(x[:6],s[:6],'o',c='r')
# plt.plot(x[:6],chisq_kept2,'o')
# plt.savefig('test.pdf')
