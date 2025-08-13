import numpy as np

class PCA:
    def __init__(self, nz_fid, pcs_path, npcs_nz):
        self.nbar = nz_fid.copy()
        self.pcs_path = pcs_path
        self.npcs_nz = npcs_nz
        self.U = np.genfromtxt(self.pcs_path)[:,:self.npcs_nz]
        self.s = self.nbar[:,1:].shape
        self.z = self.nbar[:,0]
    
    def pca(self, params_values):
        # Model: n(z) = <n>(z) + α_1*PC_1(z) + α_2*PC_2(z) + ... + α_n*PC_n(z)
        alphas = np.array([params_values.get("roman_alpha_"+str(i+1)) for i in range(self.npcs_nz)])
        correction = (alphas * self.U).sum(axis=1)

        nz_model = self.nbar[:,1:].T.flatten() + correction
        nz_model = nz_model.reshape(self.s[::-1]).T
        nz_model = np.column_stack((self.z,nz_model))
        return nz_model

# Adapted from cosmosis/samplers/fisher/* to use with CoCoA's prototype
class Fisher:
    def __init__(self,ci,nz_fid,step_size,fisher_file):
        """
        Input:
        ci: the cosmolike_{survey}_interface. e.g. roman_real
        step_size: the derivative step size
        nz_fid: the reference nz around which we want to find the derivative of the 2pt correlation function
        fisher_file: the name or path to where store the fisher matrix.
        Outputs:
        fisher matrix, jacobian matrix, and the inverse of the masked covariance matrix
        -----
        Input nz_fid has shape (Nz,1+Nt), i.e. in CoCoA .nz like-format, 
        Skip the 1st column in nz_fid (redshift). 
        Transformed nz_fid = (Nt*Nz,) =
          [ 
          nz(z1)_t1,nz(z2)_t1,nz(z3)_t1,...,nz(zNz)_t1,
          nz(z1)_t2,nz(z2)_t2,nz(z3)_t2,...,nz(zNz)_t2, 
          ..., 
          nz(z1)_Nt,nz(z2)_Nt,nz(z3)_Nt,...,nz(zNz)_Nt
          ]. 
          E.g.: For Roman, Nt=9; for DES, Nt=4 or 6; for LSST, Nt = 5.
        """
        self.ci = ci
        self.nz_fid = nz_fid.copy() # shape: (Nz,1+Nt)
        self.step_size = step_size
        self.fisher_file = fisher_file
        ################
        self.ci.set_source_sample(self.nz_fid)
        self.z = self.nz_fid[:,0]
        (self.Nz, self.Nt) = self.nz_fid[:,1:].shape
        self.tomo_bins = [(i,j) for i in range(self.Nt) for j in range(self.Nt) if j>=i] # Tomo bin combinations
        self.nparams = self.Nt * self.Nz
 
    def five_points_stencil_points(self, param_index):
        delta = np.zeros(self.nparams) # (Nt*Nz,)
        delta[param_index] = 1.0
        points_normalized=np.zeros((4,self.Nt*self.Nz))
        nz_fid_flat = self.nz_fid[:,1:].T.flatten() # shape: (Nt*Nz,) OBS: CosmoSIS perturb the normalized parameter. Here the physical parameter is pertubed directly.
        
        points = np.array([nz_fid_flat + x*delta for x in
                           [
                               +2*self.step_size, # forward far
                               +1*self.step_size, # forward near
                               -1*self.step_size, # backward near 
                               -2*self.step_size  # backward far
                           ]]) # (4, Nt*Nz)
        
        #  Normalize each bin for each stencil point
        for p in range(4):
            for t in range(self.Nt):
                start = t * self.Nz
                end = (t + 1) * self.Nz
                y = points[p, start:end]
                y /= np.trapz(y,x=self.z)
                points_normalized[p,start:end] = y

        return list(points_normalized)

    def generate_sample_points(self):
        points = []
        for p in range(self.nparams):
            points += self.five_points_stencil_points(p)
        return points

    def compute_ξpm(self):
        observable = []
        points = self.generate_sample_points() # (4*Nt*Nz,Nt*Nz)
        print('POINTS.SHAPE: ',np.array(points).shape)
        for idx,point in enumerate(points):
            print('idx: ',idx)
            point = point.reshape(self.Nt,self.Nz).T # (Nt*Nz,) -> (Nt,Nz) -> (Nz,Nt)
            point = np.column_stack((self.z,point)) # (Nz,1+Nt) CoCoA .nz like-format
            self.ci.set_source_sample(point)
            (ξ_p, ξ_m) = self.ci.xi_pm_tomo() # ξ_p (Nθ,Nt,Nt), ξ_m (Nθ,Nt,Nt)
            ξp = np.array([ξ_p[:,tbi,tbj] for (tbi,tbj) in self.tomo_bins]).flatten() # shape: Nθ*int(factorial(Nt+2-1)/(2*factorial(Nt-1)))
            ξm = np.array([ξ_m[:,tbi,tbj] for (tbi,tbj) in self.tomo_bins]).flatten() # shape: Nθ*int(factorial(Nt+2-1)/(2*factorial(Nt-1)))
            ξpm = np.hstack((ξp,ξm)) # ξpm.shape = 2 * ξp.shape
            observable.append(ξpm) 
        return observable   # (4*Nt*Nz,len(ξpm))

    def five_point_stencil_deriv(self, obs):
        obs = np.array(obs)
        deriv = (-obs[0] + 8*obs[1] - 8*obs[2] + obs[3]) / (12*self.step_size)
        return deriv

    def extract_derivatives(self, results):
        derivatives = []
        for p in range(self.nparams):
            results_p = results[4*p:4*(p+1)]
            derivative = self.five_point_stencil_deriv(results_p)
            derivatives.append(derivative)
        return np.array(derivatives)

    def central_difference(self,args):
        # truncate error ~ O(ϵ^2)
        z_idx, tomo_z_bin = args

        # ξ+_u = ξ+^{ij}(θ_k; n(z_l) + ϵ). u = up displacement
        nz_per_u = self.nz_fid.copy()
        nz_per_u[z_idx, tomo_z_bin+1] += self.step_size
        nz_per_u[:,tomo_z_bin+1] /= np.trapz(y=nz_per_u[:,tomo_z_bin+1], x=self.nz_fid[:,0])
        self.ci.set_source_sample(nz_per_u)
        ξp_per_u = self.ci.xi_pm_tomo()[0].copy() # (Nθ, Nt, Nt)
        ξm_per_u = self.ci.xi_pm_tomo()[1].copy() # (Nθ, Nt, Nt)
        
        # ξ+_d = ξ+^{ij}(θ_k; n(z_l) - ϵ). d = down displacement
        nz_per_d = self.nz_fid.copy()
        nz_per_d[z_idx, tomo_z_bin+1] -= self.step_size
        nz_per_d[nz_per_d[:, tomo_z_bin+1]<0, tomo_z_bin+1] = 0 # masks rows where column tomo_z_bin+1 is negative and use fancy indexing
        nz_per_d[:,tomo_z_bin+1] /= np.trapz(y=nz_per_d[:,tomo_z_bin+1], x=self.nz_fid[:,0])
        self.ci.set_source_sample(nz_per_d)
        ξp_per_d = self.ci.xi_pm_tomo()[0].copy() # (Nθ, Nt, Nt)
        ξm_per_d = self.ci.xi_pm_tomo()[1].copy() # (Nθ, Nt, Nt)

        # dξ+/dn = ( ξ+_p - ξ+_m ) / 2ϵ            
        dξp_dn = np.array([((ξp_per_u[:,tbi,tbj] - ξp_per_d[:,tbi,tbj]) / (2*self.step_size)) for (tbi,tbj) in self.tomo_bins]).flatten() # (len(tomo_bins),Nθ) -> (len(tomo_bins)*Nθ,). In the case of Roman (len(tomo_bins)*Nθ,)=(45*15,)=(675,): [dξp_dn^{0,0},dξp_dn^{0,1},...,dξp_dn^{0,8},dξp_dn^{1,1},...dξp_dn^{1,8},...,dξp_dn^{7,7},dξp_dn^{7,8},dξp_dn^{8,8}]
        dξm_dn = np.array([((ξm_per_u[:,tbi,tbj] - ξm_per_d[:,tbi,tbj]) / (2*self.step_size)) for (tbi,tbj) in self.tomo_bins]).flatten() # (len(tomo_bins),Nθ) -> (len(tomo_bins)*Nθ,). In the case of Roman (len(tomo_bins)*Nθ,)=(45*15,)=(675,): [dξm_dn^{0,0},dξm_dn^{0,1},...,dξm_dn^{0,8},dξm_dn^{1,1},...dξm_dn^{1,8},...,dξm_dn^{7,7},dξm_dn^{7,8},dξm_dn^{8,8}]
        
        dξpm_dn = np.hstack((dξp_dn,dξm_dn)) # (2 * len(tomo_bins)*Nθ,)
        
        return dξpm_dn

    def compute_fisher_matrix(self):

        print('COMPUTING STENCIL POINTS')
        stencil_points = self.generate_sample_points()
        np.savetxt(f'{self.fisher_file}/stencil_points_{self.step_size}.txt',stencil_points)
        print('shape: stencil_points: ', np.array(stencil_points).shape,'\n')

        print('COMPUTING XIP and XIM')
        result_ξpm = self.compute_ξpm()
        np.savetxt(f'{self.fisher_file}/ξpm_{self.step_size}.txt',result_ξpm)
        print('shape: ξpm: ', np.array(result_ξpm).shape,'\n')

        print('COMPUTING JACOBIAN MATRIX FOR XIPM')
        deriv_ξpm = self.extract_derivatives(result_ξpm)
        np.savetxt(f'{self.fisher_file}/deriv_ξpm_{self.step_size}.txt',deriv_ξpm)
        print('shape: derivative ξpm: ', np.array(deriv_ξpm).shape,'\n')
        
        # print('COMPUTING FISHER MATRIX')
        # fisher_matrix = deriv_ξpm @ self.inv_cov[:deriv_ξpm.shape[1],:deriv_ξpm.shape[1]] @ deriv_ξpm.T
        # np.savetxt(f'{self.fisher_file}/fisher_matrix_{self.step_size}.txt',fisher_matrix)
        # print('shape: fisher_matrix_dv: ', np.array(fisher_matrix).shape,'\n')
        return None

    def compute_fisher_matrix2(self):
        deriv_ξpm=[]

        jobs = [(zi, tb) for tb in range(self.Nt) for zi in range(self.Nz)]

        for job in jobs:
            print("perturbing n(zi)_bin at (zi, bin): ",job)
            deriv_ξpm.append(self.central_difference(job))
        np.savetxt(f'{self.fisher_file}/deriv_ξpm_{self.step_size}.txt',deriv_ξpm)
        return None
    
    def execute(self):
        self.compute_fisher_matrix()
        return None

    def execute2(self):
        print('COMPUTING INVERSE OF COVARIANCE MATRIX')
        self.inv_cov = self.ci.get_inv_cov_masked()
        np.savetxt(f'{self.fisher_file}/inv_cov.txt',self.inv_cov)
        print('shape: inv covariance matrix: ', np.array(self.inv_cov).shape,'\n')

        self.compute_fisher_matrix2()
        return None