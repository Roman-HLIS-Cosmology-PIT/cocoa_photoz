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

path='/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/'


# from Readme: https://docs.google.com/document/d/1iLUo9ok0hWj75-D6QhRrpidgJe4xqB0W-1sKcQu2eJU/edit?tab=t.0 
# Wide field specification (area)
# Scenario 1: (DRM-like plan) 5x140 s, dark sky, Y106+J129+H158+F184 H<24.96
# Scenario 2: (proposed wide layer) 5x91 s, bright sky, H band only, H<23.93
# Scenario 3: (proposed medium layer, if placed in an Equatorial field), 5x91 s, bright sky, Y106+J129+H158, H<24.31


# Secenarios
scenarios = [
'sc1b_d4','sc1b_d5','sc1b_d6','sc1b_d7',
'sc2b_d4','sc2b_d5','sc2b_d6','sc2b_d7',
'sc3b_d4','sc3b_d5','sc3b_d7']
scenarios_label = [
'DRM-D1','DRM-D2','DRM-D3','DRM-D4',
'W-D1','W-D2','W-D3','W-D4',
'M-D1','M-D2',        'M-D4']

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
                            sharex=True, sharey=True,
                            gridspec_kw={'wspace':0.0, 'hspace':0.0})
    axes = axes.flatten()  # flatten to easily iterate

    for ax,(sc,scl) in list(enumerate(zip(scenarios,scenarios_label))):
        nzr = f'{path}roman_nz_realizations/{sc}/nz_samples_LHC0_pointZ_1e6_Roman_{sc}.h5'
        nzr = h5py.File(nzr,'r') 
        zr = np.array(nzr['zbinsc'])
        Nz = len(zr)
        nzmax = 0
        
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

        info = {'color':colors[ax],'alpha':0.05} 
        for r in range(tot):
            [axes[ax].plot(zr,rows[r][Nz*k:Nz*(k+1)],**info,label=scl if k==0 else None) for k in range(9)]
            if r == 0:
                axes[ax].legend(loc='upper right',fontsize=17, handlelength=0, handletextpad=0)
        # Hide unused subplot
        if len(axes) > len(scenarios):
            axes[len(scenarios)].axis('off')
        # Hide tick labels on inner subplots
        n_plots=11
        ncols=4
        nrows=3
        for i in range(n_plots):
            row = i // ncols
            col = i % ncols

            # Only bottom row gets x-axis labels
            if row != nrows - 1:
                axes[i].tick_params(labelbottom=False, bottom=False)
            # Only first column gets y-axis labels
            if col != 0:
                axes[i].tick_params(labelleft=False, left=False)  
            # Only bottom row: set x-axis label
            if row == nrows - 1:
                axes[i].set_xlabel(r'$\mathrm{z}$', fontsize=20)

            # Only first column: set y-axis label
            if col == 0:
                axes[i].set_ylabel(r'$\mathrm{n(z)}$', fontsize=20)

    plt.xlim(min(zr),max(zr))
    plt.ylim(0,nzmax)
    for ax in axes.flatten():
        # Get current tick labels and set them again with larger fontsize
        ax.set_xticklabels(ax.get_xticks(), fontsize=14)
        ax.set_yticklabels(ax.get_yticks(), fontsize=14)
        # ax.yaxis.set_major_formatter(FuncFormatter(lambda val, pos: f'{int(val)}'))
        # ax.xaxis.set_major_formatter(FuncFormatter(lambda val, pos: f'{int(val)}'))
    plt.tight_layout()    
    plt.savefig(f'test.pdf')
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

def distribution_violinplot():
    plt.figure()
    
    nzs = f'{path}roman_nz_realizations/sc1b_d4/nz_samples_LHC0_pointZ_1e6_Roman_sc1b_d4.h5'
    nzs = h5py.File(nzs,'r') 
    zbins = np.array(nzs['zbinsc'])
    nzs = np.stack([nzs[f'bin0'][:100,:]], axis=1)

    plt.violinplot(nzs[:,0],positions=zbins,
                   widths=0.1, showmeans=True,
                   showmedians=True, showextrema=True)
    plt.xlabel(r'$\mathrm{z}$',fontsize=18)
    plt.ylabel(r'$\text{Distribution of } n^{tomo=1}(z)$',fontsize=18)
    plt.tight_layout()   
    plt.savefig('test.pdf')
    return None

# distribution_violinplot()