## pipeline step 1 - see pipeline_workflow.txt

## Libs
import numpy as np
import argparse
import os
import h5py

## Rule of thumb
cocoa_path = os.getcwd()
assert os.path.basename(cocoa_path) == "Cocoa", "Please run this code from ./Cocoa"

## Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--scenario', type=str, required=True, help='Roman scenario')

args = parser.parse_args()
sc = args.scenario

## Check if simulation exist
possible_scenarios = ['sc1bd4','sc1bd5','sc1bd6','sc1bd7', 'sc2bd4','sc2bd5','sc2bd6','sc2bd7', 'sc3bd4','sc3bd5','sc3bd7','sc1bd4_dz001']
if sc not in possible_scenarios: 
    raise ValueError(f"Invalid Roman scenario: '{sc}'! Should be one of: {possible_scenarios}")


file_path = f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/nz_samples_LHC0_pointZ_1e6_Roman_{sc}.h5'
nz        = h5py.File(file_path,'r') 
z         = np.array(nz['zbinsc'])
tot       = int(nz['bin0'].shape[0])
Nt        = 9                             # number of tomographic bins (presumably the same)
Nz        = int(np.min(nz['bin0'].shape)) # number of redshift bins (may change)

print(f"Number of tomographic bins {Nt}")
print(f"Number of redshift bins {Nz}")
print(f'Total simulations: {tot}')

print("nz.keys(): ",nz.keys())
print("len(z) s: ",len(np.array(nz['zbins'])))
print("len(z) center: ",len(np.array(nz['zbinsc'])))
print("len(nz['bin0'][0]): ",len(nz['bin0'][0]))

def compute_nnorm_nbar_ndiff_covn(sc=None):
    print('--------------------------------')
    print(f'PROCESSING ROMAN SCENARIO {sc}')
    print('--------------------------------')

    nzs_flatten = []

    for j in range(tot):
        print(f"Normalizing n(z) # {j}")
        row   = []
        nz_j  = np.array([nz[f'bin{i}'][j] for i in range(Nt)]).T # shape (Nz,Nt)  - Not normalized
        norms = np.trapz(y=nz_j, x=z, axis=0)                     # shape (Nt,)
        nz_j /= norms                                             # shape (Nz,Nt)  - Normalized
        nz_f  = nz_j.flatten(order="F")                           # shape (Nz*Nt,) - Flatten in column-major
        row.append(nz_f)
        nzs_flatten.append(np.hstack(row))                        # shape (tot,Nz*Nt)


    nbar_flatten  = np.mean(nzs_flatten,axis=0)                   # shape (Nz*Nt,)
    nbar_reshaped = nbar_flatten.flatten().reshape((Nt, Nz)).T
    nbar          = np.column_stack((z, nbar_reshaped))
    ndiff         = nzs_flatten - nbar_flatten
    Cn            = np.einsum('ij,ik->jk',ndiff,ndiff) / (ndiff.shape[0]-1)

    np.save(f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/normalized_nz_{sc}.npy'   , nzs_flatten)
    np.savetxt(f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/nbar_flatten_{sc}.txt' , nbar_flatten)
    np.savetxt(f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/nbar_{sc}.txt'         , nbar)
    np.save(f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/ndiff_{sc}.npy'           , ndiff)
    np.savetxt(f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/Cn_{sc}.txt'           , Cn)
    
    print()
    print("CHECKING SHAPES")
    print('n(z)s.shape:         ', np.array(nzs_flatten).shape)
    print('nbar_flatten.shape:  ', nbar_flatten.shape)
    print('nbar_reshaped.shape: ', nbar_reshaped.shape)
    print('nbar.shape:          ', nbar_flatten.shape)
    print('ndiff.shape:         ', ndiff.shape)
    print('Cn1.shape:           ', Cn.shape)

    print()
    print("FILES SAVED:")
    print('normalized n(z)s saved  @: ' ,f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/normalized_nz_{sc}.npy')
    print('mean n(z) saved         @: ' ,f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/nbar_flatten_{sc}.txt')
    print('mean n(z) saved         @: ' ,f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/nbar_{sc}.txt')
    print('difference matrix saved @: ' ,f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/ndiff_{sc}.npy')
    print('Covariance Cn saved     @: ' ,f'{cocoa_path}/cocoa_photoz/roman_nz_realizations/{sc}/Cn_{sc}.txt')
    
    return None

compute_nnorm_nbar_ndiff_covn(sc)
