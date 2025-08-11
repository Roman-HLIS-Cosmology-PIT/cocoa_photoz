import numpy as np

class PCA:
    def __init__(self, nbar_path, pcs_path, npcs_nz):
        self.nbar_path = nbar_path
        self.pcs_path = pcs_path
        self.npcs_nz = npcs_nz
        self.nbar = np.genfromtxt(self.nbar_path)
        self.U = np.genfromtxt(self.pcs_path)[:,:self.npcs_nz]
        self.s = self.nbar[:,1:].shape
        self.z = self.nbar[:,0]
    
    def pca(self, params_values):
        # Model: n(z) = <n>(z) + α_1*PC_1(z) + α_2*PC_2(z) + ... + α_n*PC_n(z)
        if self.npcs_nz > 0:
            alphas = np.array([params_values.get("roman_alpha_"+str(i+1)) for i in range(self.npcs_nz)])
            correction = (alphas * self.U).sum(axis=1)
    
            nz_model = self.nbar[:,1:].T.flatten() + correction
            nz_model = nz_model.reshape(self.s[::-1]).T
            nz_model = np.column_stack((self.z,nz_model))
            return nz_model
        else:
            return self.nbar

# Adapted from cosmosis/samplers/fisher/* to use with CoCoA's prototype
class Fisher:
    def __init__(self,ci,step_size,start_vector,fisher_file):
        """
        Input start_vector is (Nz,1+Nt) in CoCoA .nz like-format, 
        Skip the 1st column in start_vector (redshift). 
        Transformed start_vector = (Nt*Nz,) =
          [ 
          nz(z1)_t1,nz(z2)_t1,nz(z3)_t1,...,nz(zNz)_t1,
          nz(z1)_t2,nz(z2)_t2,nz(z3)_t2,...,nz(zNz)_t2, 
          ..., 
          nz(z1)_Nt,nz(z2)_Nt,nz(z3)_Nt,...,nz(zNz)_Nt
          ]. 
          E.g.: For Roman, Nt=9; for DES, Nt=4 or 6; for LSST, Nt = 5.
        """
        self.ci = ci
        self.step_size = step_size
        self.start_vector = np.genfromtxt(start_vector) # shape: (Nz,1+Nt)
        self.fisher_file = fisher_file
        ################
        self.z = self.start_vector[:,0]
        (self.Nz, self.Nt) = self.start_vector[:,1:].shape
        self.start_vector = self.start_vector[:,1:].T.flatten() # shape: (Nt*Nz,)
        self.nparams = self.Nz * self.Nt
        self.ijs = [(i,j) for i in range(self.Nt) for j in range(self.Nt) if j>=i] # Tomo bin combinations
    
    def five_points_stencil_points(self, param_index):
        delta = np.zeros(self.nparams) # (Nt*Nz,)
        delta[param_index] = 1.0
        points_normalized=np.zeros((4,self.Nt*self.Nz))
        
        points = np.array([self.start_vector + x*delta for x in
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
            ξp = np.array([ξ_p[:,ij[0],ij[1]] for ij in self.ijs]).flatten() # shape: Nθ*int(factorial(Nt+2-1)/(2*factorial(Nt-1)))
            ξm = np.array([ξ_m[:,ij[0],ij[1]] for ij in self.ijs]).flatten() # shape: Nθ*int(factorial(Nt+2-1)/(2*factorial(Nt-1)))
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
        
        print('COMPUTING INVERSE OF COVARIANCE MATRIX')
        inv_cov = self.ci.get_inv_cov_masked()
        np.savetxt(f'{self.fisher_file}/inv_cov.txt',inv_cov)
        print('shape: inv covariance matrix: ', np.array(inv_cov).shape,'\n')
        
        print('COMPUTING FISHER MATRIX')
        fisher_matrix = deriv_ξpm @ inv_cov[:deriv_ξpm.shape[1],:deriv_ξpm.shape[1]] @ deriv_ξpm.T
        np.savetxt(f'{self.fisher_file}/fisher_matrix_{self.step_size}.txt',fisher_matrix)
        print('shape: fisher_matrix_dv: ', np.array(fisher_matrix).shape,'\n')
        
        return None
    
    def execute(self):
        self.compute_fisher_matrix()
        return None