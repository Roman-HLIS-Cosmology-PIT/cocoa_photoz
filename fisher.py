import numpy as np

epsilon=0.01

class Fisher:
    def __init__(self, ci,xip_fid,n_tomo,jacob_dim1):
        self.ci = ci
        self.xip_fid = xip_fid
        self.n_tomo = n_tomo
        self.jacob_dim1 = jacob_dim1

    def compute_derivative(self,args,source_nz_local,nz_fid):
        z_idx, tomo_z_bin = args

        nz_per = source_nz_local # self.source_nz.copy()
        nz_per[z_idx, tomo_z_bin+1] += epsilon
        nz_per[:,tomo_z_bin+1] /= np.trapz(y=nz_per[:,tomo_z_bin+1], x=nz_fid[:,0])

        self.ci.set_source_sample(nz_per)
        xip_per = self.ci.xi_pm_tomo()[0].copy()
        
        dxi_dn = np.array([((xip_per[:,tbi,tbj] - self.xip_fid[:,tbi,tbj]) / epsilon) for tbi in range(self.n_tomo) for tbj in range(self.n_tomo) if tbj>=tbi]).reshape(1,self.jacob_dim1)

        return dxi_dn