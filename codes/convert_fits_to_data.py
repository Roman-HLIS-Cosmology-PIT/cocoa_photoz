"""
Extract twopt data vector, cov and nz from the fits file (cosmosis format) and 
convert to cocoa format files (cocoa/Cocoa/external_modules/data)

Applied here to FIDUCIAL_DVv5.0_lingbias27-11-24_Tz_WZ_bqr_0d01_pile3.fits
"""

###############
# Basic libs
###############
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

###############################################
## file with two-point correlation functions ##
###############################################

## CosmoSIS
fits_file = 'roman_nz_realizations/Fisher_matrix/FIDUCIAL_DVv5.0_lingbias27-11-24_Tz_WZ_bqr_0d01_pile3.fits'
## CoCoA 
pile3_path = '../projects/des_y3/data/pile3/'
data_file = 'pile3_data.txt'
nz_source_file = 'pile3_source.nz'
nz_lens_file = 'pile3_lens.nz'
cov_file = 'pile3_cov.txt'
mask_file = 'pile3_no_mask.mask'
dataset_file = 'pile3.dataset'

hdul = fits.open(fits_file)
print(hdul.info())

# source = np.array(hdul['nz_source'].data)
# lens = np.array(hdul['nz_source'].data)
# xip = np.array(hdul['xip'].data)
# gammat = np.array(hdul['gammat'].data)
# wtheta = np.array(hdul['wtheta'].data)

# print(source.dtype.names,source.shape)
# print(lens.dtype.names,lens.shape)
# print(xip.dtype.names,xip.shape)
# print(gammat.dtype.names,gammat.shape)
# print(wtheta.dtype.names,wtheta.shape)
# print(wtheta,wtheta.shape)

def modelvector():
    xip = np.array(hdul['xip'].data) # cosmic shear   
    xim = np.array(hdul['xim'].data) # cosmic shear
    gammat = np.array(hdul['gammat'].data) # galaxy galaxy lensing
    wtheta = np.array(hdul['wtheta'].data) # galaxy galaxy clustering

    xip_list = [xip[i][3] for i in range(len(xip))]
    xim_list = [xim[i][3] for i in range(len(xim))]
    gammat_list = [gammat[i][3] for i in range(len(gammat))]
    wtheta_list = [wtheta[i][3] for i in range(len(wtheta))]

    ## START ONLY XIP
    s_ntomo=4
    l_ntomo=4
    n_θ=26
    ξp_dim  = int(( s_ntomo*(s_ntomo+1) / 2 ) * n_θ)
    ξm_dim  = int((s_ntomo*(s_ntomo+1) / 2 ) * n_θ)
    gammat_dim = int((s_ntomo * l_ntomo ) * n_θ)
    w_dim = int(s_ntomo * n_θ)

    no_ξm = list(np.zeros(ξm_dim))
    no_gammat = list(np.zeros(gammat_dim))
    no_wtheta = list(np.zeros(w_dim))

    modelvector_list = enumerate(xip_list + no_ξm + no_gammat + no_wtheta)
    np.savetxt(data_file,list(modelvector_list), fmt='%d %.18e')
    ## END ONLY XIP

    # modelvector_list = enumerate(xip_list + xim_list + gammat_list + wtheta_list)
    # np.savetxt(data_file,list(modelvector_list), fmt='%d %.18e')

    return None
# modelvector()

def nz_source_lens():
    nz_source = np.vstack(hdul['nz_source'].data)[1:, [1, 3, 4, 5, 6]].copy() # source
    nz_lens = np.vstack(hdul['nz_lens'].data)[1:, [1, 3, 4, 5, 6]].copy() # lens
    np.savetxt(nz_source_file,nz_source, fmt='%.2f %.18e %.18e %.18e %.18e')
    np.savetxt(nz_lens_file,nz_lens, fmt='%.2f %.18e %.18e %.18e %.18e')
    return None
# nz_source_lens()


def covariance():
    # covmat = np.array(hdul['COVMAT'].data)
    # dim = covmat.shape[0]
    # covmat = covmat.reshape(int(covmat.shape[0]**2),-1)
    # bin_comb = [[i,j] for i in range(dim) for j in range(dim)]
    # covmat = [bin_comb[i]+list(covmat[i]) for i in range(dim**2)]
    # np.savetxt(cov_file,covmat,fmt='%d %d %.18e')
    
    ## START ONLY XIP
    covmat=np.array(hdul['COVMAT'].data)
    s_ntomo=4
    l_ntomo=4
    n_θ=26
    ξp_dim  = int(( s_ntomo*(s_ntomo+1) / 2 ) * n_θ)
    ξm_dim  = int((s_ntomo*(s_ntomo+1) / 2 ) * n_θ)
    gammat_dim = int((s_ntomo * l_ntomo ) * n_θ)
    w_dim = int(s_ntomo * n_θ)
    ndata_ = ξp_dim+ξm_dim+gammat_dim+w_dim ## dimension that cosmolike wants: [critical] this->ndata_ = 1040
    #### TODO START ####
    # new_covmat = np.zeros((ndata_, ndata_)) ## RuntimeError: inv(): matrix is singular
    # new_covmat = np.ones((ndata_, ndata_))*1e-20 ## [critical] IP::set_inv_cov: masked cov not positive definite
    #### TODO END ####
    new_covmat[:ξp_dim, :ξp_dim] = covmat[:ξp_dim, :ξp_dim]
    new_covmat = new_covmat.reshape(int(new_covmat.shape[0]**2),-1)
    bin_comb = [[i,j] for i in range(ndata_) for j in range(ndata_)]
    new_covmat = [bin_comb[i]+list(new_covmat[i]) for i in range(ndata_**2)]
    np.savetxt(cov_file,new_covmat,fmt='%d %d %.18e')
    ## END ONLY XIP
    return None
covariance()


def no_mask():
    # dv_size = np.array(hdul['COVMAT'].data).shape[0]
    # ones = [(i,1.0) for i in range(dv_size)]
    # np.savetxt(mask_file,ones,fmt='%d %.1f')

    ## START ONLY XIP
    s_ntomo=4
    l_ntomo=4
    n_θ=26
    ξp_dim  = int(( s_ntomo*(s_ntomo+1) / 2 ) * n_θ)
    ξm_dim  = int((s_ntomo*(s_ntomo+1) / 2 ) * n_θ)
    gammat_dim = int((s_ntomo * l_ntomo ) * n_θ)
    w_dim = int(s_ntomo * n_θ)
    ndata_ = ξp_dim+ξm_dim+gammat_dim+w_dim ## cosmolike wants this dimension: [critical] this->ndata_ = 1040
    ones = [(i,1.0) for i in range(ndata_)]
    np.savetxt(mask_file,ones,fmt='%d %.1f')
    ## END ONLY XIP    
    return None
# no_mask()

def dataset():
    xip = np.array(hdul['xip'].data)
    theta_min_arcmin = min(xip['ANGLEMIN'])
    theta_max_arcmin = max(xip['ANGLEMAX'])

    template=f"""data_file = {data_file}
cov_file = {cov_file}
mask_file = {mask_file}
nz_lens_file = {nz_lens_file}
nz_source_file = {nz_source_file}
lens_ntomo = 4
source_ntomo = 4
n_theta = 26
theta_min_arcmin = {theta_min_arcmin}
theta_max_arcmin = {theta_max_arcmin}
baryon_pca_file = pca.txt"""

    with open(dataset_file,'w') as f:
        f.write(template)
    return None    
# dataset()

def execute_all():
    modelvector()
    print("Model vector - DONE")
    nz_source_lens()
    print("Source and lens - DONE")
    covariance()
    print("Covariance - DONE")
    no_mask()
    print("Mask - DONE")
    dataset()
    print("Dataset - DONE")
    return None
# execute_all()