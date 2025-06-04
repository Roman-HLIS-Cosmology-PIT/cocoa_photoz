import numpy as np
import matplotlib.pyplot as plt

path = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/"

eps_values1 = [1e-20,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,2e-1,3e-1,4e-1,5e-1,6e-1,7e-1]

row = 0
col = 0
temp = []

def figure1():
    plt.figure()
    for i,eps in enumerate(eps_values):
        x = np.genfromtxt(path+f"test_central_difference_eps{eps_values1[i]}.txt")[row,col]
        print(x)
        temp.append(x)

    plt.plot(eps_values,abs(np.array(temp)))
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\epsilon$",fontsize=15)
    plt.ylabel(r"$\partial\xi_{+}^{ij}(\theta_l)~/~\partial n_k$",fontsize=15)
    plt.tight_layout()
    plt.savefig("./figures/epsilon_stability.pdf")

figure1()


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

# plt.figure()
# deriv2 = np.genfromtxt(path+f"test_central_difference_eps{eps_values[-1]}.txt")
# im2 = plt.imshow(deriv2,cmap="seismic")
# plt.colorbar(im2,orientation='vertical')
# plt.tight_layout()
# plt.savefig("deriv2.pdf")