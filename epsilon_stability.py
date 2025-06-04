import numpy as np
import matplotlib.pyplot as plt

path = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/"

# eps_values1 = [1e-20,...,1e-1,2e-1,3e-1,4e-1,5e-1,6e-1,7e-1]
log_a = -14
eps_values1 = list(np.round(np.logspace(log_a,-1,-log_a),decimals=20))+[2e-1,3e-1,4e-1,5e-1,6e-1,7e-1]

row = 0
col = 0

dxipdn_forward = []
dxipdn_central = []

def figure1():
    plt.figure()
    for i,eps in enumerate(eps_values1):
        x1 = np.genfromtxt(path+f"test_forward_difference/test_forward_difference_eps{eps_values1[i]}.txt")[row,col]
        x2 = np.genfromtxt(path+f"test_central_difference/test_central_difference_eps{eps_values1[i]}.txt")[row,col]
        print(x1,x2)
        dxipdn_forward.append(x1)
        dxipdn_central.append(x2)

    plt.plot(eps_values1,abs(np.array(dxipdn_central)),c='k',ls='-',label='central difference')
    plt.plot(eps_values1,abs(np.array(dxipdn_forward)),c='C0',ls='--',label='forward difference')
    plt.xscale("log")
    plt.yscale("log")
    plt.title('LSST Y1 Source')
    plt.xlabel(r"$\epsilon$",fontsize=15)
    plt.ylabel(r"$|~\partial\xi_{+}^{ij}(\theta_a)~/~\partial n^k(z_b)~|$",fontsize=15)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("./epsilon_stability.pdf")

# figure1()


eps_values = [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,2e-1,3e-1,4e-1,5e-1,6e-1,7e-1]
def figure2():
    plt.figure()

    fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(10, 5))
    ii = 0
    for i in range(4):
        for j in range(4):        
            print(ii)
            deriv = np.genfromtxt(path+f"test_central_difference_eps{eps_values[ii]}.txt")
            im = axes[i,j].imshow(deriv,cmap="seismic",aspect='auto')
            plt.colorbar(im,orientation='vertical')
            axes[i,j].set_title(r"$\epsilon=$"+f"{eps_values[ii]}")
            if i == 3 and j == 0:
                axes[i,j].set_xlabel(r"$\theta$",fontsize=15)
                axes[i,j].set_ylabel(r"$n(z)$",fontsize=15)
            ii+=1

    plt.tight_layout()
    plt.savefig("./figures/deriv.pdf")

# plt.figure2()