import numpy as np

covdim = 2700

def no_mask():
    no_mask = list(enumerate(np.ones(covdim)))
    np.savetxt('no_mask.mask',no_mask, fmt='%d %.1f')
    return None