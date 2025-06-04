import numpy as np
import matplotlib.pyplot as plt


cov_lsst_y1 = np.genfromtxt("/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/lsst_y1/data/lsst_y1_cov")
cov_des_y3 = np.genfromtxt("/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/des_y3/data/des_y3_cov_unblinded_final.txt")
cov_roman_real = np.genfromtxt("/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/roman_real/data/cov_roman_real")

print(cov_des_y3.shape,cov_des_y3.shape[0]**.5)
print(cov_lsst_y1.shape,cov_lsst_y1.shape[0]**.5)
print(cov_roman_real.shape,cov_roman_real.shape[0]**.5)



# dim = (x.shape[0])**0.5

# t = x.shape[0]
# print(np.sqrt(t))

# cov_g00 = x[0:26*26,8]
# dim_cov_g00 = int(len(x[0:26*26,8])**(0.5))

# cov_g00 = cov_g00.reshape(dim_cov_g00,dim_cov_g00)

# im = plt.imshow(cov_g00,cmap="seismic")
# plt.colorbar(im,orientation='vertical')
# plt.tight_layout()
# plt.savefig("./lsst_y1_cov.pdf")
