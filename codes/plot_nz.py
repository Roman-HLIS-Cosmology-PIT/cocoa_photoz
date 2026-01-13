#!/usr/bin/env python
#SBATCH --job-name=plot_nz
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=7-00:00:00

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib

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

path='/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa/cocoa_photoz/'


# from Readme: https://docs.google.com/document/d/1iLUo9ok0hWj75-D6QhRrpidgJe4xqB0W-1sKcQu2eJU/edit?tab=t.0 
# Wide field specification (area)
# Scenario 1: (DRM-like plan) 5x140 s, dark sky, Y106+J129+H158+F184 H<24.96
# Scenario 2: (proposed wide layer) 5x91 s, bright sky, H band only, H<23.93
# Scenario 3: (proposed medium layer, if placed in an Equatorial field), 5x91 s, bright sky, Y106+J129+H158, H<24.31

# Secenarios
scenarios = [
'sc1bd4','sc1bd5','sc1bd6','sc1bd7',
'sc2bd4','sc2bd5','sc2bd6','sc2bd7',
'sc3bd4','sc3bd5','sc3bd7']
scenarios_label = [
'M1-D1','M1-D2','M1-D3','M1-D4',
'M2-D1','M2-D2','M2-D3','M2-D4',
'M3-D1','M3-D2',        'M3-D4']

# A thousand random integers between 1 and 1M
np.random.seed(42)
sims_indexes = np.random.randint(1, 10**6 + 1, size=100)
tot = len(sims_indexes)

# Colors
colors = [
    "#1b5f6f",  # dark teal
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#999999",  # grey
    "#117733",  # dark green
    "#882255",  # wine
]

# Plot the thousands nzs
def plot_scenarios(sc='sc1b_d4'):

    print(f'Total simulations: {tot}')

    # Create a figure with subplots (3 rows)
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(14, 8),
                            sharex=False, sharey=False,
                            gridspec_kw={'wspace':0.0, 'hspace':0.0})
    axes = axes.flatten()  # flatten to easily iterate
 
    for ax,(sc,scl) in list(enumerate(zip(scenarios,scenarios_label))):
        nzr = f'{path}roman_nz_realizations/{sc}/nz_samples_LHC0_pointZ_1e6_Roman_{sc}.h5'
        nzr = h5py.File(nzr,'r') 
        zr = np.array(nzr['zbinsc'])
        Nz = len(zr)
        nzmax = 0
        nbar = np.genfromtxt(f'{path}roman_nz_realizations/{sc}/nbar_{sc}.nz')
        
        print(f'Nz: {Nz}, Scenario: {sc}, Axes: {ax}')

        rows = []
        for j in sims_indexes:
            # print(f"Processing j={j}")
            row = []
            for i in range(9):
                normalized = nzr[f'bin{i}'][j] / np.trapz(y=nzr[f'bin{i}'][j], x=zr)
                row.append(normalized)
            rows.append(np.hstack(row))
        
        nzmax_temp = np.max(np.array(rows))
        if nzmax_temp > nzmax:
            nzmax = nzmax_temp

        info = {"color":colors[ax],"alpha":0.03,"lw":1} 
        temp=0
        # The nz realizations
        for r in range(tot):
            if ax==10: temp=1
            [axes[ax+temp].plot(zr,rows[r][Nz*k:Nz*(k+1)],**info,label=scl if k==0 else None) for k in range(9)]
            if r == 0:
                axes[ax+temp].legend(loc='upper right',fontsize=17, handlelength=0, handletextpad=0)
        # The mean nz
        for i in range(1,9+1):
            axes[ax+temp].plot(zr,nbar[:,i],color=colors[ax],lw=1.5,ls="--")        

    # Hide the axis of scenario M3-D3 
    axes[10].axis('off')

    for i in [6,8,9,10,11]:
        axes[i].set_xlabel(r'$z$', fontsize=20)
        axes[i].set_xticks([0,0.5,1,1.5,2,2.5])
    for i in [0,4,8]:
        axes[i].set_ylabel(r'$n(z)$', fontsize=20)
        axes[i].set_yticks([0,1,2,3,4])
    for i in [1,2,3,5,6,7,9,10,11]:
        axes[i].set_yticks([])
    for i in range(11):
        if i != 8:
            axes[i].set_xlim(np.min(zr),np.max(zr))
    for i in range(11):
        axes[i].set_ylim(0,nzmax)
    # plt.tight_layout()    
    fig_name = "roman_scenarios.pdf"
    plt.savefig(f"{fig_name}")
    print(f"Figure {fig_name}")
    # n = np.vstack(rows)
    # print('Saving n')
    # print('n.shape: ',n.shape)
    # # np.savetxt('n_roman_sc1bd4.txt',n)
    # np.save('n_roman_sc1bd4.npy',n)

    # nbar = np.mean(n,axis=0)
    # print('Saving nbar')
    # print('nbar.shape: ',nbar.shape)
    # # np.savetxt('nbar_roman_sc1bd4.txt',nbar)
    # np.save('nbar_roman_sc1bd4.npy',nbar)

    # ndiff = n - nbar
    # print('Saving ndiff')
    # print('ndiff.shape: ',ndiff.shape)
    # # np.savetxt('ndiff_roman_sc1bd4.txt',ndiff)
    # np.save('ndiff_roman_sc1bd4.npy',ndiff)
    return None
plot_scenarios()

def plot_ndiff():
    nzr = f'{path}roman_nz_realizations/sc1b_d4/nz_samples_LHC0_pointZ_1e6_Roman_sc1b_d4.h5'
    nzr = h5py.File(nzr,'r') 
    zr = np.array(nzr['zbinsc'])
    ndiff = np.load('../ndiff_roman_sc1bd4.npy')
    nbar = np.load('../nbar_roman_sc1bd4.npy')
    nrows = ndiff.shape[0] # 1M
    ncols = ndiff.shape[1] # 414
    Nz = len(zr) # 46

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(5, 4),
                            sharex=True, sharey=False,
                            gridspec_kw={'hspace':0.0, 'height_ratios':[1,0.3]})

    info = {'color':colors[0],'alpha':0.05} 
    [axes[0].plot(zr,nbar[Nz*k:Nz*(k+1)],color='gray',alpha=0.6) for k in range(9)]
    for idx in sims_indexes:
        [axes[1].plot(zr,ndiff[idx,Nz*k:Nz*(k+1)],**info) for k in range(9)]
    axes[0].set_xlim(np.min(zr),np.max(zr)) 
    axes[0].set_ylim(np.min(nbar),np.max(nbar)) 
    axes[0].set_ylabel(r'$\bar{\mathbf{n}}_\text{ref}$',fontsize=18)
    axes[1].set_xlabel(r'$\mathrm{z}$',fontsize=18)
    axes[1].set_ylabel(r'$\mathbf{n}_\text{ref}-\bar{\mathbf{n}}_\text{ref}$',fontsize=18)
    plt.tight_layout()   
    plt.savefig('test.pdf')    
    return None

# plot_ndiff()

def distribution_violinplot(sc='sc1bd4'):
    #####
    fig, axes = plt.subplots(nrows=3, ncols=9, figsize=(14, 5),
                        sharex=False, sharey=True,
                        gridspec_kw={'wspace':0.0, 'hspace':0.0})

    Nsim=5000
    indices=[8,13,17,19,22,25,29,32,37]
    slices = [(0,4),(4,8),(8,11)]

    temp=0
    for row,s in enumerate(slices): 
        left,right = s[0],s[1]
        for i,sc in enumerate(scenarios[left:right]):
            print(sc)
            nzsf = f'{path}roman_nz_realizations/{sc}/nz_samples_LHC0_pointZ_1e6_Roman_{sc}.h5'
            nbar = np.genfromtxt(f'{path}roman_nz_realizations/{sc}/nbar_{sc}.nz')
            nzs = h5py.File(nzsf,'r') 
            zbins = np.array(nzs['zbinsc'])
            for t,ind in enumerate(indices):
                arr = nzs[f'bin{t}'][:Nsim, :]
                arr_norm = arr / np.trapz(arr, x=zbins, axis=1)[:, None]
                nzs_t = np.stack([arr_norm], axis=1)
                zbar = np.trapz(y=nbar[:,0][:,None]*nbar[:,1:],x=nbar[:,0],axis=0)/\
                    np.trapz(y=nbar[:,1:],x=nbar[:,0],axis=0)

                vp=axes[row,t].violinplot(dataset=nzs_t[:,0][:,ind],positions=[zbins[ind]],
                            widths=0.1, showmeans=False,
                            showmedians=False, showextrema=False)

                for body in vp['bodies']:
                    body.set_facecolor('none')
                    body.set_edgecolor(colors[i+temp])     
                    body.set_linewidth(1)
                    body.set_alpha(1)     
        temp+=4
    #####    

    # plt.scatter(np.full_like(nzs[:,0][:,8], zbins[8]),nzs[:,0][:,8], marker='o', linestyle='None',s=4)
    axes[1,0].set_ylabel(r'$\text{Distribution}$',fontsize=18)
    axes[2,4].set_xlabel(r'$\mathrm{arg\,min}_{z_i} \, |z_i - \bar{z}|$', fontsize=18)
    plt.tight_layout()   
    plt.savefig('distribution_violinplot.pdf')
    return None

# distribution_violinplot()

def nz_Delta_LensingEfficiencyKernel():

    import camb.constants
    from scipy.interpolate import interp1d
    from scipy.integrate import quad

    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
    results = camb.get_background(pars)
    nbar = np.loadtxt("cocoa_photoz/roman_nz_realizations/sc1bd4/nbar_sc1bd4.nz")

    z_tab = nbar[:,0]
    const_factor = 3 * 0.5 * pars.H0**2 * pars.omegam / (camb.constants.c*1e-3)**2

    def Wi(z_tab,nzi_tab):
        """
        z_tab: tabulated redshifts
        nzi_tab: tabulated redshift distribution at bin i
        """
        chi = results.comoving_radial_distance(z_tab)

        dchidz = np.gradient(chi,z_tab)

        int_Wi_1 = nzi_tab * dchidz
        int_Wi_2 = nzi_tab * dchidz / chi
        int_Wi_1_interp = interp1d(z_tab, int_Wi_1)
        int_Wi_2_interp = interp1d(z_tab, int_Wi_2)
        pre_fac = chi*(1+z_tab)
        pre_fac_interp = interp1d(z_tab, pre_fac, kind='cubic', bounds_error=False, fill_value='extrapolate')

        Wi_1 = lambda z: quad(int_Wi_1_interp, z, z_tab[-1])[0]
        Wi_2 = lambda z: quad(int_Wi_2_interp, z, z_tab[-1])[0]
        Wi_temp = lambda z: const_factor*pre_fac_interp(z)*(Wi_1(z) - results.comoving_radial_distance(z)*Wi_2(z))
        return np.array([Wi_temp(zi) for zi in z_tab],dtype=np.float64)

    for i in range(1,10):
        plt.plot(z_tab,Wi(z_tab,nbar[:,i]))

    sc = 'sc1bd4'
    cocoa_path = '/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa' # Change for your path 
    nz = h5py.File(f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/nz_samples_LHC0_pointZ_1e6_Roman_{sc}.h5','r')
    nz_mean = np.genfromtxt(f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/nbar_sc1bd4.nz')
    z  = np.array(nz['zbinsc'])
    Nt = 9
    Nz = 46
    nz_tomos=[]
    Nsim = 100
    for r in range(Nsim):
        nz_tomos.append([np.column_stack(nz[f'bin{i}'][r]/np.trapz(y=nz[f'bin{i}'][r],x=z))[0] for i in range(9)])

    nz_tomos=np.array(nz_tomos) # shape (Nsim,Nt,Nz)
    nz_tomos=np.transpose(nz_tomos, (0, 2, 1)) # shape (Nsim,Nz,Nt)
    print(nz_tomos.shape)
    print(nz_tomos[0,:,0].shape)

    # plt.figure()
    # fig = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=False,figsize=(8,3),
    #                         gridspec_kw={'wspace':0.2, 'hspace':0.0,'height_ratios': [2, 1]})

    fig = plt.figure(figsize=(8, 6))

    gs = fig.add_gridspec(2, 2, wspace=0.4, hspace=0.3)

    # Left column: two separate subplots
    axes1 = fig.add_subplot(gs[0, :])  # top-left
    axes2 = fig.add_subplot(gs[1, 0])  # bottom-left
    axes3 = fig.add_subplot(gs[1, 1])  # all rows, column 1

    for t in range(Nt):
        print('bin:',t)
        for s in range(Nsim):
            print(s,end=',')
            axes1.plot(z,nz_tomos[s,:,t],c='#1b5f6f',alpha=.03)
            axes3.plot(z,Wi(z,nz_tomos[s,:,t]),c='#1b5f6f',alpha=.03)

    for i in range(9): 
        axes1.plot(z,nz_mean[:,i+1],c='#1b5f6f',ls='--')

    for r in range(Nsim):
        for i in range(9):
            axes2.plot(z,nz_tomos[r,:,i]-nz_mean[:,i+1],c='gray',alpha=.03)


    axes1.set_ylabel(r'$\bar{\mathbf{n}}$', labelpad=1)
    axes2.set_ylabel(r'${\mathbf{\Delta}}$', labelpad=1)
    axes3.set_ylabel(r'$W_i$', labelpad=1)
    axes1.set_xlabel(r'$z$')
    axes2.set_xlabel(r'$z$')
    axes3.set_xlabel(r'$z$')
    axes1.set_yticks([1,2,3]);
    axes1.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0]);
    axes2.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0]);
    axes3.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0]);
    axes1.set_xlim(np.min(z),np.max(z));
    axes2.set_xlim(np.min(z),np.max(z));
    axes2.set_xlim(np.min(z),np.max(z));
    axes1.set_ylim(0,np.max(nz_mean))
    axes2.set_ylim(-1,+1)
    axes2.set_yticks([-1,-0.5, 0,0.5,1])
    # fig.align_ylabels(axes1)
    # axes[1,0].grid(False)
    plt.tight_layout()
    plt.savefig('redshift_delta_sc1bd4_v2.pdf');    
    return None

# nz_Delta_LensingEfficiencyKernel()