import numpy as np
import matplotlib.pyplot as plt
from getdist import plots,loadMCSamples,MCSamples

burnin=0.5
settings = {"ignore_rows": burnin}
cp = '../results/chains' # common path
mcmc_control = {'c_1':{'path':f'{cp}/roman_sc1bd4_g/',
                       'prefix':'PCA_MODEL',
                       'info':'8 decimal in nz. Status: (canceled)'},
                'c_2':{'path':f'{cp}/roman_sc1bd4_g_18decimal/',
                       'prefix':'PCA_MODEL',
                       'info':'18 decimal in nz. Status: (canceled)'},
                'c_3':{'path':f'{cp}/roman_sc1bd4_test_dv/',
                       'prefix':'PCA_MODEL_modelvector.py',
                       'info':'dv generated with modelvector.py. Status: (canceled)'},
                'c_4':{'path':f'{cp}/same_params_than_cosmocov/',
                       'prefix':'PCA_MODEL',
                       'info':'dv generated with modelvector.py. Status: (running)'},
                'c_5':{'path':f'{cp}/diff_params_than_cosmocov/',
                       'prefix':'PCA_MODEL',
                       'info':'dv generated with modelvector.py. Status: (running)'},
                'c_6':{'path':f'{cp}/roman_real/',
                       'prefix':'ROMAN_REAL',
                       'info':'example1.modelvector: 3x2 (original)'},
                'c_7':{'path':f'{cp}/diff_params_than_cosmocov_gauss_alphas_prior/',
                       'prefix':'PCA_MODEL',
                       'info':'Same as c_5, but with Gaussian prior on alphas.'},
                       }

control = 'c_7'
path = mcmc_control[control]['path']
px = mcmc_control[control]['prefix']
info = mcmc_control[control]['info']
print(info)

def roman_pca_mcmc():
#     no_pc = loadMCSamples(path+f"{px}_MCMC0", settings=settings)
    one_pc = loadMCSamples(path+f"{px}_MCMC1", settings=settings)
    two_pc = loadMCSamples(path+f"{px}_MCMC2", settings=settings)
    three_pc = loadMCSamples(path+f"{px}_MCMC3", settings=settings)

    g = plots.get_subplot_plotter()
#     g.plots_1d([no_pc,one_pc,two_pc,three_pc],["omegam","sigma8","w","wa"],
#             markers={"omegam":0.3,"sigma8":0.8120,"w":-1,"wa":0}, # planck best-fit 1807.06209
#             nx=4,legend_ncol=5,legend_labels=['No PC', '1 PC', '2 PC', '3 PC'],
#             colors=['#1b5f6f',"#E69F00","#D55E00","#56B4E9"],
#             ls=['-','--','-.',':'],lws=[3,3,3,3])
    g.plots_1d([one_pc,two_pc,three_pc],["omegam","sigma8","w","wa"],
            markers={"omegam":0.3,"sigma8":0.8120,"w":-1,"wa":0}, # planck best-fit 1807.06209
            nx=4,legend_ncol=5,legend_labels=['1 PC', '2 PC', '3 PC'],
            colors=["#E69F00","#D55E00","#56B4E9"],
            ls=['--','-.',':'],lws=[3,3,3])
    g.export("1d_marg.pdf")

    g = plots.get_subplot_plotter()
    g.plots_1d([one_pc],["roman_alpha_1"],
            markers={"roman_alpha_1":0.0},
            colors=["#E69F00"],
            ls=['-'],lws=[3])
    g.export("1d_marg_alpha_1pc.pdf")

    g = plots.get_subplot_plotter()
    g.plots_1d([two_pc],["roman_alpha_1","roman_alpha_2"],
            markers={"roman_alpha_1":0.0,"roman_alpha_2":0.0},
            colors=["#D55E00"],nx=2,
            ls=['-'],lws=[3])
    g.export("1d_marg_alpha_2pc.pdf")

    g = plots.get_subplot_plotter()
    g.plots_1d([three_pc],["roman_alpha_1","roman_alpha_2","roman_alpha_3"],
            markers={"roman_alpha_1":0.0,"roman_alpha_2":0.0,"roman_alpha_3":0.0},
            colors=["#56B4E9"],nx=3,
            ls=['-'],lws=[3])
    g.export("1d_marg_alpha_3pc.pdf")
    return None
roman_pca_mcmc()

def roman_real_mcmc():
    mcmc0 = loadMCSamples(path+f"{px}_MCMC0", settings=settings)
    mcmc1 = loadMCSamples(path+f"{px}_MCMC1", settings=settings)
    mcmc2 = loadMCSamples(path+f"{px}_MCMC2", settings=settings)

    g = plots.get_subplot_plotter()
    g.plots_1d([mcmc0,mcmc1,mcmc2],["omegam","sigma8","w","wa"],
            markers={"omegam":0.3,"sigma8":0.8120,"w":-1,"wa":0}, # planck best-fit 1807.06209
            nx=4,legend_ncol=5,
            legend_labels=['ROMAN_REAL_MCMC0', 'ROMAN_REAL_MCMC1', 'ROMAN_REAL_MCMC2'],
            colors=['#1b5f6f',"#E69F00","#D55E00"],
            ls=['-','--','-.'],lws=[3,3,3])
    g.export("1d_marg.pdf")

    g = plots.get_subplot_plotter()
    g.plots_1d([mcmc0,mcmc1,mcmc2],["roman_DZ_S1","roman_DZ_S2","roman_DZ_S3"],
            markers={"roman_DZ_S1":0.0,"roman_DZ_S2":0.0,"roman_DZ_S3":0.0},
            colors=['#1b5f6f',"#E69F00","#D55E00"],nx=3,
            ls=['-','--','-.'],lws=[3,3,3])
    g.export("1d_marg_alpha_3pc.pdf")
    return None

def plot_RMinus1():
    plt.figure()
    Rminus1_stop = 0.015
    max_val = 0
    progress_0 = f'{px}_MCMC0.progress'
    progress_1 = f'{px}_MCMC1.progress'
    progress_2 = f'{px}_MCMC2.progress'
    progress_3 = f'{px}_MCMC3.progress'

#     progress = [progress_0,progress_1,progress_2,progress_3]
#     colors = ['#1b5f6f',"#E69F00","#D55E00","#56B4E9"]
#     legend_labels=['No PC', '1 PC', '2 PC', '3 PC']
    
    progress = [progress_1,progress_2,progress_3]
    colors = ["#E69F00","#D55E00","#56B4E9"]
    legend_labels=['1 PC', '2 PC', '3 PC']

    # progress = [progress_0,progress_1,progress_2]
    # colors = ['#1b5f6f',"#E69F00","#D55E00"]
    # legend_labels=['ROMAN_REAL_MCMC0', 'ROMAN_REAL_MCMC1', 'ROMAN_REAL_MCMC2']

    for p,c,l in zip(progress,colors,legend_labels):
        Rminus1 = np.genfromtxt(path+p,usecols=(3))
        Rminus1 = Rminus1[~np.isnan(Rminus1)]
        vals = list(range(len(Rminus1)))
        if  np.max(vals) > max_val:
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