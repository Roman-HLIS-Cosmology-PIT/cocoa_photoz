import numpy as np
import matplotlib.pyplot as plt
import h5py

path='/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/'
nz = f'{path}/roman_nz_realizations/sc1b_d4/SVSN/nz_samples__LHC0_pointZ_1e6_Roman_sc1b_d4.h5'
nz = h5py.File(nz,'r') 
z = np.array(nz['zbinsc'])

def plot_modes():
    U = np.load(f'{path}/U_source.npz')
    nEig = U['U'].shape[0]  # number of eigenvectors
    ntomo = U['U'].shape[1] # number of tomographic bins
    zbins = U['U'].shape[2] # number of redhisft bins

    c=['r','orange','darkgreen','darkblue']

    fig,ax = plt.subplots(3,3,figsize=(10,10))

    grid=[(i,j) for i in range(3) for j in range(3)]
    grid = list(enumerate(grid))

    for t, (row,col) in grid:
        for eig in range(4):
            mode = U['U'][eig,t,:]
            config = {'c':c[eig],
                    'label':[f'{eig}'],
                    'ls':'--' if eig==len(range(4))-1 else None}
            ax[row,col].plot(z,mode,**config)
    ax[0,0].legend(loc='best')
    plt.savefig('U_source.pdf')
    return None
# plot_modes()

def plot_chisqkept():
    chisq_kept = np.genfromtxt(f'{path}/chisq_kept.txt')
    Ms = np.arange(len(chisq_kept))
    plt.plot(Ms,chisq_kept,'o')
    plt.savefig('chisq_kept.pdf')
    return None
plot_chisqkept()