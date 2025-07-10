import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

dvpile3 = './FIDUCIAL_DVv5.0_lingbias27-11-24_Tz_WZ_bqr_0d01_pile3.fits'
hdul = fits.open(dvpile3)
covmat = np.array(hdul['COVMAT'].data)
covmat_reduced = covmat[0:260,0:260]

plt.figure()
# im_full = plt.imshow(np.linalg.inv(covmat),cmap='seismic',vmin=-1,vmax=+1)
im_full = plt.imshow(covmat,cmap='seismic',vmin=-1,vmax=+1)
plt.colorbar(im_full)
plt.savefig('./inv_covmat_full_pile3.pdf')

# covmat_3x2 = np.array(hdul['COVMAT'].data)

# covdim = covmat_3x2.shape[0]


# print(covmat_3x2.shape)

# pp_norm = np.zeros((covdim,covdim))
# for i in range(covdim):
#     for j in range(covdim):
#         pp_norm[i][j] = covmat_3x2[i][j] / np.sqrt(covmat_3x2[i][i]*covmat_3x2[j][j])

# invcov = np.linalg.inv(pp_norm)

# plt.imshow(invcov,cmap="seismic",vmin=-1,vmax=+1)
# plt.colorbar()

# plt.savefig("./des_pile3_cov.pdf")