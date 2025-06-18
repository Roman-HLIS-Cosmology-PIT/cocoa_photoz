import numpy as np
import matplotlib.pyplot as plt
import h5py

path='/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/'
nz = f'{path}/roman_nz_realizations/sc1b_d4/SVSN/nz_samples__LHC0_pointZ_1e6_Roman_sc1b_d4.h5'
nz = h5py.File(nz,'r') 
z = np.array(nz['zbinsc'])

def plot_modes(thr):
    plt.figure()
    U = np.load(f'{path}/U_source_{thr}.npz')
    nEig = U['U'].shape[0]  # number of eigenvectors
    ntomo = U['U'].shape[1] # number of tomographic bins
    zbins = U['U'].shape[2] # number of redhisft bins

    # c=['r','orange','darkgreen','darkblue']

    fig,ax = plt.subplots(3,3,figsize=(13,13),sharex=True,sharey=False)

    grid=[(i,j) for i in range(3) for j in range(3)]
    grid = list(enumerate(grid))

    for t, (row,col) in grid:
        for eig in range(nEig):
            mode = U['U'][eig,t,:]
            config = {
                    # 'c':c[eig],
                    'label':[f'{eig}'],
                    'ls':'--' if eig==len(range(nEig))-1 else None}
            ax[row,col].plot(z,mode,**config)
    ax[0,0].legend(loc='upper right')
    ax[2,0].set_xlabel('z',fontsize=12)
    ax[2,1].set_xlabel('z',fontsize=12)
    ax[2,2].set_xlabel('z',fontsize=12)
    ax[0,0].set_ylabel('Modes',fontsize=12)
    ax[1,0].set_ylabel('Modes',fontsize=12)
    ax[2,0].set_ylabel('Modes',fontsize=12)
    ax[0,1].set_title(r'$\langle\chi^2_{th}\rangle=$'+f'{thr}', fontsize=16)
    plt.savefig(f'U_source_{thr}.pdf')
    return None

thr = [0.001,0.005,0.01]
for t in thr:
    plot_modes(t)

def plot_chisqkept():
    plt.figure()
    thr = [0.001,0.005,0.01]
    chisq_kept0 = np.genfromtxt(f'{path}/chisq_kept_{thr[0]}.txt')
    chisq_kept1 = np.genfromtxt(f'{path}/chisq_kept_{thr[1]}.txt')
    chisq_kept2 = np.genfromtxt(f'{path}/chisq_kept_{thr[2]}.txt')
    Ms0 = np.arange(len(chisq_kept0))
    Ms1 = np.arange(len(chisq_kept1))
    Ms2 = np.arange(len(chisq_kept2))
    c = ['purple','C0','orange']
    plt.plot(Ms0,chisq_kept0,'o',c=c[0])
    plt.plot(Ms1-0.2,chisq_kept1,'s',c=c[1])
    plt.plot(Ms2+0.2,chisq_kept2,'^',c=c[2])
    minv = min(min(Ms0),min(Ms1),min(Ms2))
    maxv = max(max(Ms0),max(Ms1),max(Ms2))
    info = lambda i: {'xmin':minv-.3,'xmax':maxv,
            'color':f'{c[i]}',
            'ls':':',
            'label':r'$\langle\chi^2_{th}\rangle=$'+f'{thr[i]}'}
    plt.hlines(thr[0],**info(0))
    plt.hlines(thr[1],**info(1))
    plt.hlines(thr[2],**info(2))
    plt.xticks(range(minv,maxv + 1, 1))
    plt.xlim(minv-.3,maxv)
    plt.xlabel('M',fontsize=13)
    plt.ylabel(r'$\langle\chi^2_{th}\rangle$',fontsize=13)
    plt.legend(loc='best')
    plt.savefig(f'chisq_kept.pdf')
    return None
# plot_chisqkept()