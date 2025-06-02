import numpy as np
from math import factorial # DHFS MOD

epsilon=0.01

class Fisher:
    def __init__(self,ci,nz_fid,n_tomo,n_theta):
        self.ci = ci
        self.nz_fid = nz_fid
        self.n_tomo = n_tomo
        self.n_theta = n_theta
        self.jacob_dim1 = self.n_theta * int(factorial(self.n_tomo+2-1)/(2*factorial(self.n_tomo-1)))
        self.ci.set_source_sample(self.nz_fid)
        self.xip_fid = self.ci.xi_pm_tomo()[0].copy()

    def forward_difference(self,args):
        # error ~ ϵ
        z_idx, tomo_z_bin = args

        nz_per = self.nz_fid.copy()
        nz_per[z_idx, tomo_z_bin+1] += epsilon
        nz_per[:,tomo_z_bin+1] /= np.trapz(y=nz_per[:,tomo_z_bin+1], x=self.nz_fid[:,0])

        self.ci.set_source_sample(nz_per)
        xip_per = self.ci.xi_pm_tomo()[0].copy()
        
        dxi_dn = np.array([((xip_per[:,tbi,tbj] - self.xip_fid[:,tbi,tbj]) / epsilon) for tbi in range(self.n_tomo) for tbj in range(self.n_tomo) if tbj>=tbi]).reshape(1,self.jacob_dim1)

        return dxi_dn
    
    def central_difference(self,args):
        # error ~ ϵ^2
        z_idx, tomo_z_bin = args

        nz_per_p = self.nz_fid.copy()
        nz_per_p[z_idx, tomo_z_bin+1] += epsilon
        nz_per_p[:,tomo_z_bin+1] /= np.trapz(y=nz_per_p[:,tomo_z_bin+1], x=self.nz_fid[:,0])
        self.ci.set_source_sample(nz_per_p)
        xip_per_p = self.ci.xi_pm_tomo()[0].copy()
        
        nz_per_m = self.nz_fid.copy()
        nz_per_m[z_idx, tomo_z_bin+1] -= epsilon
        nz_per_m[nz_per_m[:, tomo_z_bin+1]<0, tomo_z_bin+1] = 0 # maks rows where column tomo_z_bin+1 is negative and use fancy indexing
        nz_per_m[:,tomo_z_bin+1] /= np.trapz(y=nz_per_m[:,tomo_z_bin+1], x=self.nz_fid[:,0])
        self.ci.set_source_sample(nz_per_m)
        xip_per_m = self.ci.xi_pm_tomo()[0].copy()
        
        dxi_dn = np.array([((xip_per_p[:,tbi,tbj] - xip_per_m[:,tbi,tbj]) / (2*epsilon)) for tbi in range(self.n_tomo) for tbj in range(self.n_tomo) if tbj>=tbi]).reshape(1,self.jacob_dim1)

        return dxi_dn