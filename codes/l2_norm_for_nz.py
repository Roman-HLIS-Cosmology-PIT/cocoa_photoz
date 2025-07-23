# Code to find the nz that differ the most from the its mean nz
import numpy as np
import h5py

path='/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa'
sc='sc1bd4'
params_values={}

nz = h5py.File(f'{path}/cocoa_photoz/roman_nz_realizations/{sc}/nz_samples_LHC0_pointZ_1e6_Roman_{sc}.h5','r')
nbar_path = f'{path}/projects/roman_real/data/sc1bd4/nbar_sc1bd4.nz'

z = nz['zbinsc'][:]
nbar = np.genfromtxt(nbar_path)[:,1:] # (46,9)
Nt = 9
Nz = 46
tot_sims = nz['bin0'].shape[0]

# 10 different ways we can classify the 'biggest distance^2' in the 46-D vectorial space
biggest_norm2 = -1
biggest_norm2_index = -1
biggest_norm2_bin = np.ones(Nt)*-1
biggest_norm2_bin_indexes = np.ones(Nt,dtype=int)*-1

for sim in range(tot_sims):
    nz_normalized = []
    print(f'PROCESSING REALIZATION # {sim}')
    for i in range(Nt):
        nzi = nz[f'bin{i}'][sim,:]
        nzi /= np.trapz(y=nzi,x=z)
        nz_normalized.append(nzi)
    nz_normalized=np.array(nz_normalized).T # (46,9)

    ndiff2 = (nz_normalized - nbar)**2
    norm2_per_bin = np.sum(ndiff2,axis=0) # (9,) sum along rows - compress columns \equiv column-wise sum
    norm2_total = np.sum(norm2_per_bin)
    if norm2_total > biggest_norm2:
        biggest_norm2 = norm2_total
        biggest_norm2_index = sim

    for j in range(Nt):
        if norm2_per_bin[j] > biggest_norm2_bin[j]:
            biggest_norm2_bin[j] = norm2_per_bin[j]
            biggest_norm2_bin_indexes[j] = sim

biggest_norms2_final = list(biggest_norm2_bin_indexes)
biggest_norms2_final.insert(0,biggest_norm2_index)

np.savetxt(f'{path}/cocoa_photoz/codes/biggest_norms2_final.txt',biggest_norms2_final, fmt='%d')

print(biggest_norm2_index)
print(biggest_norm2_bin_indexes)
print(biggest_norms2_final)
print('END OF COMPUTATION')