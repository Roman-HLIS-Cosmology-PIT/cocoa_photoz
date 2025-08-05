import numpy as np
import matplotlib.pyplot as plt
import h5py

nzr0 = np.genfromtxt('../../projects/roman_real/data/sc1bd4/nz_sc1bd4r0.nz')
nbar = np.genfromtxt('../../projects/roman_real/data/sc1bd4/nbar_sc1bd4.nz')
nbar2 = np.genfromtxt('../n_nbar_ndiff_covndiff/nbar_sc1b_d4.nz')
nz = h5py.File('../roman_nz_realizations/sc1b_d4/nz_samples_LHC0_pointZ_1e6_Roman_sc1b_d4.h5','r')

# print(nzr0.shape)          # (46,10)
# print(nbar.shape)          # (46,10)
# print(nz['bin0'][0].shape) # (46,)

# print(nbar.shape)
# print(nbar2.shape)

for i in range(9):
    # print(nz[f'bin{i}'][0]/np.trapz(y=nz[f'bin{i}'][0],x=nzr0[:,0])==nzr0[:,i+1])
    print(nbar[:,i]==nbar2[:,i])
    print("")