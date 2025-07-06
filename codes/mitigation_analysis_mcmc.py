import numpy as np
import matplotlib.pyplot as plt
from getdist import plots,loadMCSamples,MCSamples

burnin=0.5
settings = {"ignore_rows": burnin}
# path = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/results/chains/roman_sc1bd4_g/" # 8 decimal in nz (stoped)
# path = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/results/chains/roman_sc1bd4_g_18decimal/" # 18 decimal in nz (stoped)
# path = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/results/chains/roman_sc1bd4_test_dv/" # dv generated with modelvector.py
# path = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/results/chains/same_params_than_cosmocov/" # dv generated with modelvector.py
path = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/results/chains/diff_params_than_cosmocov/" # dv generated with modelvector.py
# path1 = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/results/chains/roman_sc1bd4_g/ROMAN_REAL_MCMC0_example1.modelvector" # example1.modelvector: 3x2 (original)
# path2 = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/results/chains/roman_sc1bd4_g/ROMAN_REAL_MCMC0_print_datavector_example1.modelvector" # cosmic shear only from print data vector
# path3 = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/results/chains/roman_sc1bd4_g/ROMAN_REAL_MCMC0_print_datavector_example1_3x2pt.modelvector" # 3x2 from print data vector
px = 'PCA_MODEL'
# px = 'PCA_MODEL_modelvector.py'
# px = 'PCA_MODEL_printdv_3x2'

# no_pc = loadMCSamples(path+f"{px}_MCMC0", settings=settings)
# one_pc = loadMCSamples(path+f"{px}_MCMC1", settings=settings)
# two_pc = loadMCSamples(path+f"{px}_MCMC2", settings=settings)
# three_pc = loadMCSamples(path+f"{px}_MCMC3", settings=settings)


# g = plots.get_subplot_plotter()
# g.plots_1d([no_pc,one_pc,two_pc,three_pc],["omegam","sigma8","w","wa"],
#            markers={"omegam":0.3,"sigma8":0.8120,"w":-1,"wa":0}, # planck best-fit 1807.06209
#            nx=4,legend_ncol=5,legend_labels=['No PC', '1 PC', '2 PC', '3 PC'],
#            colors=['#1b5f6f',"#E69F00","#D55E00","#56B4E9"],
#            ls=['-','--','-.',':'],lws=[3,3,3,3])
# g.export("1d_marg.pdf")

def plot_RMinus1():
    plt.figure()
    Rminus1_stop = 0.015
    max_val = 0
    progress_0 = 'PCA_MODEL_MCMC0.progress'
    progress_1 = 'PCA_MODEL_MCMC1.progress'
    progress_2 = 'PCA_MODEL_MCMC2.progress'
    progress_3 = 'PCA_MODEL_MCMC3.progress'
    progress = [progress_0,progress_1,progress_2,progress_3]
    colors = ['#1b5f6f',"#E69F00","#D55E00","#56B4E9"]
    legend_labels=['No PC', '1 PC', '2 PC', '3 PC']

    for p,c,l in zip(progress,colors,legend_labels):
        Rminus1 = np.genfromtxt(path+p,usecols=(3))
        Rminus1 = Rminus1[~np.isnan(Rminus1)]
        vals = list(range(len(Rminus1)))
        if np.max(vals) > max_val:
            max_val = np.max(vals)
        plt.plot(vals,Rminus1,label=l,marker='o',markersize=3,color=c)

    plt.hlines(Rminus1_stop,np.min(vals),max_val,color='gray',ls='--',label=f'Rminus_stop: {Rminus1_stop}',lw=3)
    plt.yscale('log')
    plt.legend(loc='best',fontsize=13)
    plt.ylabel('R-1',fontsize=20)
    plt.xlabel('Points',fontsize=20)
    plt.xlim(np.min(vals),max_val)
    plt.ylim(Rminus1_stop-0.1,np.max(Rminus1)+5)
    plt.savefig('Rminus1.pdf')
plot_RMinus1()