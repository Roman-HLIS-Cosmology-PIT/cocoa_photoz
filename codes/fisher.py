import numpy as np
from math import factorial

# from cobaya.likelihoods.base_classes import DataSetLikelihood
# class TEST_FISHER(DataSetLikelihood):

# epsilon=0.01

class Fisher:
    def __init__(self,ci,nz_fid,n_tomo,n_theta,epsilon):
        self.ci = ci
        self.nz_fid = nz_fid
        self.n_tomo = n_tomo
        self.n_theta = n_theta
        self.epsilon = epsilon
        self.jacob_dim1 = self.n_theta * int(factorial(self.n_tomo+2-1)/(2*factorial(self.n_tomo-1)))
        self.ci.set_source_sample(self.nz_fid)
        self.ξp_fid = self.ci.xi_pm_tomo()[0].copy()
        self.ξm_fid = self.ci.xi_pm_tomo()[1].copy()
        self.inv_cov = self.ci.get_inv_cov_masked()

    def forward_difference(self,args):
        # truncate error ~ O(ϵ)
        z_idx, tomo_z_bin = args

        # ξ+_p = ξ+^{ij}(θ_k; n(z_l) + ϵ)
        nz_per = self.nz_fid.copy()
        nz_per[z_idx, tomo_z_bin+1] += self.epsilon
        nz_per[:,tomo_z_bin+1] /= np.trapz(y=nz_per[:,tomo_z_bin+1], x=self.nz_fid[:,0])

        # ξ+ = ξ+^{ij}(θ_k; n(z_l))
        self.ci.set_source_sample(nz_per)
        ξp_per = self.ci.xi_pm_tomo()[0].copy()
        
        # dξ_+/dn = ( ξ+_p - ξ+ ) / ϵ
        dξp_dn = np.array([((ξp_per[:,tbi,tbj] - self.ξp_fid[:,tbi,tbj]) / self.epsilon) for tbi in range(self.n_tomo) for tbj in range(self.n_tomo) if tbj>=tbi]).reshape(1,self.jacob_dim1)

        return dξp_dn
    
    def central_difference(self,args):
        # truncate error ~ O(ϵ^2)
        z_idx, tomo_z_bin = args

        # ξ+_u = ξ+^{ij}(θ_k; n(z_l) + ϵ). u = up displacement
        nz_per_u = self.nz_fid.copy()
        nz_per_u[z_idx, tomo_z_bin+1] += self.epsilon
        nz_per_u[:,tomo_z_bin+1] /= np.trapz(y=nz_per_u[:,tomo_z_bin+1], x=self.nz_fid[:,0])
        self.ci.set_source_sample(nz_per_u)
        ξp_per_u = self.ci.xi_pm_tomo()[0].copy()
        
        # ξ+_d = ξ+^{ij}(θ_k; n(z_l) - ϵ). d = down displacement
        nz_per_d = self.nz_fid.copy()
        nz_per_d[z_idx, tomo_z_bin+1] -= self.epsilon
        nz_per_d[nz_per_d[:, tomo_z_bin+1]<0, tomo_z_bin+1] = 0 # masks rows where column tomo_z_bin+1 is negative and use fancy indexing
        nz_per_d[:,tomo_z_bin+1] /= np.trapz(y=nz_per_d[:,tomo_z_bin+1], x=self.nz_fid[:,0])
        self.ci.set_source_sample(nz_per_d)
        ξp_per_d = self.ci.xi_pm_tomo()[0].copy()
        
        # dξ+/dn = ( ξ+_p - ξ+_m ) / 2ϵ
        dξp_dn = np.array([((ξp_per_u[:,tbi,tbj] - ξp_per_d[:,tbi,tbj]) / (2*self.epsilon)) for tbi in range(self.n_tomo) for tbj in range(self.n_tomo) if tbj>=tbi]).reshape(1,self.jacob_dim1)

        return dξp_dn
    
    def fisher_matrix(self,args):
        dxipdn = self.central_difference(args)
        inv_cov_xip = self.inv_cov[0:390,0:390]
        fisher_mat = dxipdn @ inv_cov_xip @ dxipdn.T
        return fisher_mat
