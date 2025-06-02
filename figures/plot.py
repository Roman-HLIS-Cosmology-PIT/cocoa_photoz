import numpy as np
import matplotlib.pyplot as plt

# path_jacob = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz"
path_jacob = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa"
derivs_cocoa = np.loadtxt(path_jacob+"/test.txt")

im = plt.imshow(derivs_cocoa,cmap="seismic")
plt.colorbar(im,orientation='vertical')
plt.tight_layout()
plt.savefig(path_jacob+"/jacobian_test.pdf")