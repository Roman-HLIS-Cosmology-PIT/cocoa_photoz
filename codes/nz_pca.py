import numpy as np

class PCA:
    def __init__(self,nz_fid,survey,npc):
        self.nz_fid = nz_fid
        self.survey = survey
        self.npc = npc
    
    def pca(self,params_values,U):
        if self.npc > 0:
            #TODO: optimization: if U is in the shape of n, reshape n will not be necessary
            alphas = np.array([params_values.get(self.survey+"_alpha_"+str(i+1)) for i in range(self.npc)])
            correction = (alphas * U).sum(axis=1)
    
            s = self.nz_fid[:,1:].shape
            z = self.nz_fid[:,0]

            nz_model = self.nz_fid[:,1:].T.flatten() + correction
            nz_model = nz_model.reshape(s[::-1]).T
            nz_model = np.column_stack((z,nz_model))
            return nz_model
        else:
            return self.nz_fid