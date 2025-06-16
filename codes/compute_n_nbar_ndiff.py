
import numpy as np
import h5py
import matplotlib.pyplot as plt

path='/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/'

def compute_n_nbar_ndiff():
    nzr = f'{path}/roman_nz_realizations/sc1b_d4/SVSN/nz_samples__LHC0_pointZ_1e6_Roman_sc1b_d4.h5'
    nzr = h5py.File(nzr,'r') 
    zr = np.array(nzr['zbinsc'])

    tot = int(nzr['bin0'].shape[0])

    print(f'Total simulations: {tot}')

    rows = []
    for j in range(tot):
        print(f"Processing j={j}")
        row = []
        for i in range(9):
            normalized = nzr[f'bin{i}'][j] / np.trapz(y=nzr[f'bin{i}'][j], x=zr)
            row.append(normalized)
        rows.append(np.hstack(row))

    n = np.vstack(rows)
    print('Saving n')
    print('n.shape: ',n.shape)
    # np.savetxt('n_roman_sc1bd4.txt',n)
    np.save('n_roman_sc1bd4.npy',n)

    nbar = np.mean(n,axis=0)
    print('Saving nbar')
    print('nbar.shape: ',nbar.shape)
    # np.savetxt('nbar_roman_sc1bd4.txt',nbar)
    np.save('nbar_roman_sc1bd4.npy',nbar)

    ndiff = n - nbar
    print('Saving ndiff')
    print('ndiff.shape: ',ndiff.shape)
    # np.savetxt('ndiff_roman_sc1bd4.txt',ndiff)
    np.save('ndiff_roman_sc1bd4.npy',ndiff)
    return None
# compute_n_nbar_ndiff()

def compute_Cn():
    print('Computing Cn - Covariance matrix of ndiff = n - nbar')
    ndiff = np.load(f'{path}/ndiff_roman_sc1bd4.npy')
    Cn = np.einsum('ij,ik->jk',ndiff,ndiff) / ndiff.shape[0]
    np.save('Cn_roman_sc1bd4.npy',Cn)
    np.savetxt('Cn_roman_sc1bd4.txt',Cn)
    return None
# compute_Cn()

# nbar  = np.genfromtxt(f'{path}/nbar_roman_sc1bd4.txt')
# print('nbar.shape: ',nbar.shape)
# ndiff = np.genfromtxt(f'{path}/ndiff_roman_sc1bd4.txt')
# print('ndiff.shape: ',ndiff.shape)

# z and nzs
# nzd = f'{path}/roman_nz_realizations/Fisher_matrix/Tz_realizations_WZ_bq_pile3_0d01.npy'
# nzd = np.load(nzd) ## shape = (10095, 4, 300) = (Ns, Nt, Nz)
# N = (3-0)/0.01
# zd = np.linspace(0,3,int(N-1))

# nzd = nzd[:,:,1:]
# nzd = nzd.reshape(np.shape(nzd)[0], 4*299)
# nbar = np.mean(nzd, axis=0)

# yr = np.array([nzr[f'bin0'][0]/np.trapz(y=nzr[f'bin0'][0],x=zr),
#                nzr[f'bin1'][0]/np.trapz(y=nzr[f'bin1'][0],x=zr),
#                nzr[f'bin2'][0]/np.trapz(y=nzr[f'bin2'][0],x=zr),
#                nzr[f'bin3'][0]/np.trapz(y=nzr[f'bin3'][0],x=zr)
#                ])

# yr2 = np.array([np.hstack(nzr[f'bin{i}'][0]) for i in range(4)])

# yr3 = np.hstack([nzr[f'bin{i}'][0]/np.trapz(y=nzr[f'bin{i}'][0],x=zr) for i in range(9)])
# yr4 = np.hstack([nzr[f'bin{i}'][1]/np.trapz(y=nzr[f'bin{i}'][1],x=zr) for i in range(9)])
# p1  = np.vstack((yr3,yr4))

# p2 = np.vstack([np.hstack([nzr[f'bin{i}'][j]/np.trapz(y=nzr[f'bin{i}'][j],x=zr) for i in range(9)]) for j in range(tot)])

# print(yr.shape,yr2.shape,yr3.shape,9*46,p1.shape,p2.shape)
# print(np.mean(p1,axis=0)==np.mean(p2,axis=0))
# print(np.mean(p2,axis=0)==np.mean(p3,axis=0))

# plt.figure()
# for i in range(9):
#     yr = nzr[f'bin{i}'][0]
#     ir = np.trapz(y=yr,x=zr)
#     yr = yr/ir
#     ir = np.trapz(y=yr,x=zr)
#     plt.plot(zr,yr)
#     print(f'∫n^{i}dz = {ir}')
# plt.savefig('nzr.pdf')

# nbar = np.mean(n, axis=0)
