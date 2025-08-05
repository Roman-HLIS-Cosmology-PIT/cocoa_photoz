import numpy as np
from itertools import chain

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
    def __init__(self,start_vector,ci):
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
        self.start_vector = np.genfromtxt(start_vector) # shape=(Nz,1+Nt)
        self.z = self.start_vector[:,0]
        self.Nz = len(self.z)
        self.Nt = self.start_vector[:,1].shape[1]
        self.start_vector = list(chain.from_iterable(start_vector[:,1].T)) # shape=(Nt*Nz,)
        self.nparams = len(self.start_vector) # = self.Nz * self.Nt
        self.ci = ci
    
    def five_points_stencil_points(self, param_index):
        delta = np.zeros(self.nparams) # Nt*Nz
        delta[param_index] = 1.0
        points = [self.start_vector + x*delta for x in
                  [
                      +2*self.step_size, # forward far
                      +1*self.step_size, # forward near
                      -1*self.step_size, # backward near 
                      -2*self.step_size  # backward far
                  ]
                 ] # (4, Nt*Nz)
        return points

    def generate_sample_points(self):
        points = []
        for index in range(self.nparams):
            points += self.five_points_stencil_points(index) # List concatenation
        return points # (4*Nt*Nz, Nt*Nz)

    def compute_obs(self):
        points = self.generate_sample_points()

        for point in points:
            point=point.copy()
            self.ci.set_source_sample(point)
            (ξp, ξm) = self.ci.xi_pm_tomo()
            (ntheta, ntomo, ntomo2) = ξp.shape
            ξp_list += ξp[:,tbi,tbj] list(chain.from_iterable(ξp))
            ξp = [ξp[:,tbi,tbj] for tbi in range(ntomo) for tbj in range(ntomo) if tbj>=tbi].reshape(1,self.jacob_dim1)
            ξp_list=[]
            ξp_list += ξp[:,]
            ξm = self.ci.xi_pm_tomo()[1].copy()

    def five_point_stencil_deriv(self, obs):
        deriv = (-obs[0] + 8*obs[1] - 8*obs[2] + obs[3]) / (12*self.step_size)
        return deriv

    def extract_derivatives(self, results):
        derivatives = []
        #Now get out the results that correspond to each dimension
        for p in range(self.nparams):
            results_p = results[4*p:4*(p+1)]
            derivative = self.five_point_stencil_deriv(results_p, p)
            derivatives.append(derivative)
        return np.array(derivatives)
    
    def compute_derivatives(self):
            derivatives = []
            points = self.generate_sample_points()
            print("Calculating derivatives using {} total models".format(len(points)))
            if self.pool is None:
                results = list(map(self.compute_vector, points))
            else:
                results = self.pool.map(self.compute_vector, points)

            _, inv_cov = self.compute_vector(points[0], cov=True)

            derivatives = self.extract_derivatives(results)

            return derivatives, inv_cov # ORIGINAL
    
    def compute_fisher_matrix(self):
        derivatives, inv_cov = self.compute_derivatives()
        fisher_matrix = np.einsum("il,lk,jk->ij", derivatives, inv_cov, derivatives)
        return fisher_matrix
    