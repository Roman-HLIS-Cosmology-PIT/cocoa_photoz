
import numpy as np
import h5py
import matplotlib.pyplot as plt
import argparse

path='/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/'
p = f'{path}n_nbar_ndiff_covndiff/'

scenarios = [
'sc1b_d4','sc1b_d5','sc1b_d6','sc1b_d7',
'sc2b_d4','sc2b_d5','sc2b_d6','sc2b_d7',
'sc3b_d4','sc3b_d5','sc3b_d7']

def compute_n_nbar_ndiff(sc=''):
    print('--------------------------------')
    print(f'PROCESSING ROMAN SCENARIO {sc}')
    print('--------------------------------')

    file = f'{path}/roman_nz_realizations/{sc}/nz_samples_LHC0_pointZ_1e6_Roman_{sc}.h5'
    nz = h5py.File(file,'r') 
    z = np.array(nz['zbinsc'])

    tot = int(nz['bin0'].shape[0])

    print(f'Total simulations: {tot}')

    rows = []
    for j in range(tot):
        print(f"Normalizing n(z) # {j}")
        row = []
        for i in range(9):
            normalized = nz[f'bin{i}'][j] / np.trapz(y=nz[f'bin{i}'][j], x=z)
            row.append(normalized)
        rows.append(np.hstack(row))

    n = np.vstack(rows)
    print('Saving n(z) normalized')
    print('n(z)_normalized.shape: ',n.shape)
    np.save(f'{p}nz_normalized_{sc}.npy',n)

    nbar = np.mean(n,axis=0)
    print('Saving mean n(z)')
    print('nbar.shape: ',nbar.shape)
    np.save(f'{p}nbar_{sc}.npy',nbar)
    reshaped_nbar = nbar.flatten().reshape((9, 46)).T
    final_nbar = np.column_stack((z, reshaped_nbar))
    np.savetxt(f'{p}nbar_{sc}.nz', final_nbar, fmt=['%.3f'] + ['%.8e'] * 9)  # Change fmt if needed for different precision

    ndiff = n - nbar
    print('Saving difference matrix')
    print('ndiff.shape: ',ndiff.shape)
    np.save(f'{p}ndiff_{sc}.npy',ndiff)
    return None
# compute_n_nbar_ndiff()

def compute_Cn(sc='sc1bd4'):
    print('Computing Cn - Covariance matrix of the difference matrix: ndiff = n - nbar')
    ndiff = np.load(f'{p}ndiff_{sc}.npy')
    Cn = np.einsum('ij,ik->jk',ndiff,ndiff) / (ndiff.shape[0]-1)
    np.save(f'{p}cov_ndiff_{sc}.npy',Cn)
    np.savetxt(f'{p}cov_ndiff_{sc}.txt',Cn)
    return None
# compute_Cn()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sc", type=str, required=True)
    args = parser.parse_args()
    compute_n_nbar_ndiff(sc=args.sc)
    compute_Cn(sc=args.sc)