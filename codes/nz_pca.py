import numpy as np

class PCA:
    def __init__(self, nbar_path, pcs_path, npcs_nz):
        self.nbar_path = nbar_path
        self.pcs_path = pcs_path
        self.npcs_nz = npcs_nz
        self.nbar = np.genfromtxt(self.nbar_path)
        self.U = np.genfromtxt(self.pcs_path)[:,:self.npcs_nz]
    
    def pca(self, params_values):
        if self.npcs_nz > 0:
            alphas = np.array([params_values.get("roman_alpha_"+str(i+1)) for i in range(self.npcs_nz)])
            correction = (alphas * self.U).sum(axis=1)
    
            s = self.nbar[:,1:].shape
            z = self.nbar[:,0]

            nz_model = self.nbar[:,1:].T.flatten() + correction
            nz_model = nz_model.reshape(s[::-1]).T
            nz_model = np.column_stack((z,nz_model))
            return nz_model
        else:
            return None