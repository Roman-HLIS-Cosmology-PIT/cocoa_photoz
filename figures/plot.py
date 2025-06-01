import numpy as np
import matplotlib.pyplot as plt

path_jacob = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz"
derivs_cocoa = np.loadtxt(path_jacob+"/results_jacobian/jacobian.txt")

im = plt.imshow(derivs_cocoa,cmap="seismic")
plt.colorbar(im,orientation='vertical')
plt.tight_layout()
plt.savefig(path_jacob+"/figures/jacobian.pdf")