import numpy as np
import matplotlib.pyplot as plt
import h5py


surname = 'fisher2_xip_no_scale_cut' ## fisher matrix

#### TO DO #####
### DHFS - REMOVE ALLLLLLL MAGIC NUMBERS! PUT NAMES!!!
#### TO DO #####

U = np.load('./U_source.npz')

print(U['U'].shape)
print(U['U'][0])
print(U['U'][0].shape)

### DES CASE ###
N = (3-0)/0.01
# z = np.linspace(0,3,int(N)-1)
z = np.linspace(0,3,int(N))

print(len(z))

fig, axs = plt.subplots(nrows=2,ncols=4,sharey=False,sharex=False,figsize=(12,7))


# U_fisher_pile3_3x2 = np.load('/gpfs/projects/MirandaGroup/Diogo/ROMAN-NZ-PROJECT/PCA/photoz_uz/NZ_COMPRESSION_OUTS/U_source_with_fisher_DES.npz')
# U_fisher_pile3_xip = np.load('/gpfs/projects/MirandaGroup/Diogo/ROMAN-NZ-PROJECT/PCA/photoz_uz/NZ_COMPRESSION_OUTS/U_source__%s.npz' % surname)

# Us = [U_fisher_pile3_3x2, U_fisher_pile3_xip]

Us = [np.load('./U_source.npz')]

for k in np.arange(1):
    U = Us[k]
    for j in np.arange(4):
        for i in np.arange(5):
            if k == 0:
                axs[k,j].plot(z,U['U'][i][j],label=rf'$\mathrm{{PC}}_{i}$' if j==0 else None)
            elif k == 1:    
                axs[k,j].plot(z,U['U'][i][j],label=rf'$\mathrm{{PC}}_{i}$' if j==0 else None)
            axs[k,j].set_xlabel('z',fontsize=15)
            if k == 0:
                axs[k,j].set_title(f'Tomo. Bin: {j}' if j==0 else f'{j}')
        axs[k,j].set_xlim(0,3)
        axs[k,j].set_ylim(-0.2,0.2)

    if k == 0:
        axs[k,0].legend(loc='upper right')
        axs[k,0].set_ylabel('CosmoSIS\nDV pile3 w/ 3x2 pt',fontsize=15)
    elif k == 1:    
        axs[k,0].legend(loc='upper right')
        axs[k,0].set_ylabel('CosmoSIS\nDV pile3 w/ '+r'$\xi_+$',fontsize=15)

plt.tight_layout()
plt.savefig('./test.pdf')