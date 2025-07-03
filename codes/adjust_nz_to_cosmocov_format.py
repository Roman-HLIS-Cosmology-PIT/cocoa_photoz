# given a flatten nz or column nz but without the redshift at the first column, we transform the nz in the format to run cosmocov

import numpy as np
import matplotlib.pyplot as plt
import h5py

nz = np.genfromtxt('../nbar_roman_sc1bd4.txt')
nz = nz.reshape(46,9,order='F')
file = h5py.File('../roman_nz_realizations/sc1b_d4/SVSN/nz_samples__LHC0_pointZ_1e6_Roman_sc1b_d4.h5','r')
z = file["zbinsc"][:]
zbinsmin = file["zbins"][:46] # z bin min (cosmocov) 
zbinsmin_reshaped = zbinsmin.reshape(-1, 1)



test = np.hstack((zbinsmin_reshaped,nz))

np.savetxt("mean_roman_nzmin_normalized_sc1bd4_40.nz", test) # z min + nzs (cosmocov)


x = np.genfromtxt('mean_roman_nzmin_normalized_sc1bd4_40.nz')
xz = x[:,0]


for i in range(1,9):
    plt.plot(z,x[:,i])
    print(i,np.trapz(y=x[:,i], x=z))
plt.savefig('test.pdf')