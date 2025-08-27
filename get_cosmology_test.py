import numpy as np
import h5py
import matplotlib.pyplot as plt
import os

# 'sc1bd4' # Already exist

scs = [
    'sc1bd5', 'sc1bd6', 'sc1bd7',
    'sc2bd4', 'sc2bd5', 'sc2bd6', 'sc2bd7',
    'sc3bd4', 'sc3bd5', 'sc3bd7'
    ]

for sc in scs:
    # sc = scs[0]
    print(sc)
    path = f'roman_nz_realizations/{sc}/nz_samples_LHC0_pointZ_1e6_Roman_{sc}.h5'
    nz = h5py.File(path)

    z = np.array(nz['zbinsc'])[:]
    Nt = 9
    Nz = 46
    nzr0 = np.zeros((Nz,1+Nt)) # 1 column for the redshift

    nzr0[:,0] = z

    for i in range(1,10):
        nzr0[:,i]=nz[f'bin{i-1}'][0]/np.trapz(y=nz[f'bin{i-1}'][0],x=z)

    os.mkdir(f'../projects/roman_real/data/{sc}')
    np.savetxt(f"../projects/roman_real/data/{sc}/nz_{sc}r0.nz",nzr0)


    # for i in range(1,10):
    #     plt.plot(z,x[:,i],'C0',ls='-')
    #     plt.plot(z,y[:,i],'darkorange',ls='--')
    #     print(np.trapz(y=x[:,i],x=z))
    #     print(np.trapz(y=y[:,i],x=z))
    # plt.savefig('test.pdf')