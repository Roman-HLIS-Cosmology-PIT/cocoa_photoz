import numpy as np
import matplotlib.pyplot as plt
from getdist import plots,loadMCSamples,MCSamples

burnin=0.5
settings = {"ignore_rows": burnin}
path = '../results/chains/roman_pca/MCMC'
# path = '../results/chains/edge_roman_pca/MCMC'
# path = '../results/datachallenge_1/MCMC'

colors = [
    '#1f77b4',  # blue
    '#ff7f0e',  # orange
    '#2ca02c',  # green
    '#d62728',  # red
    '#9467bd',  # purple
    '#8c564b',  # brown
    '#e377c2',  # pink
    '#7f7f7f',  # gray
    '#bcbd22',  # olive
    '#17becf',  # cyan
    '#aec7e8'   # light blue
]


start,end = 15,23
# idxes = [34, 35]#+list(range(start,end+1))
# idxes = [0, 12, 2, 11] #, 16, 17, 18, 19, 20, 21, 22, 23]
idxes = [0, 12] #, 16, 17, 18, 19, 20, 21, 22, 23]
print(idxes)

def roman_pca_mcmc():
    chains=[]
    for i in idxes:
        chains.append(loadMCSamples(f"{path}{i}", settings=settings)) 

    # g = plots.get_subplot_plotter()

    # g.triangle_plot(chains,
    #                 ["omegam","sigma8",
    #                 #  "roman_DZ_S1","roman_DZ_S2","roman_DZ_S3","roman_DZ_S4","roman_DZ_S5","roman_DZ_S6","roman_DZ_S7","roman_DZ_S8","roman_DZ_S9"
    #                  ],
    #                  markers={"omegam": 0.3156, "omegab": 0.0492, "sigma8": 0.812229, "H0": 67.32,}, # planck best-fit 1807.06209
    #                 #  markers={"omegam": 0.3, "omegab": 0.04, "sigma8": 0.8277, "H0": 67.32,}, # datachallenge 1
    #                  legend_ncol=1,
    #                 #  legend_labels=[],
    #                 #  line_args=[],
    #                  contour_colors=colors,
    #                  filled=True)
    # g.export("triangle_plot.pdf")
    
    # g = plots.get_subplot_plotter()
    # g.plots_1d(chains,["omegam"],
    #         markers={"omegam": 0.3156},
    #         # colors=["#D55E00"],
    #         nx=1,
    #         # ls=['-'],
    #         # lws=[3]
    #         )
    g = plots.get_subplot_plotter(width_inch=4)
    print(g,type(g))
    g.plot_1d(chains,"omegam")
    g.add_x_marker(0.3156,lw=1.5)
    g.settings.legend_fontsize=9
    # g.add_legend(["No PC + Shift + "+r"$\bar{n}$","1 PC + "+r"$\bar{n}$","10 PC + "+r"$\bar{n}$", "No PC"],legend_loc='center right',legend_ncol=1);
    g.add_legend(["No PC + Shift + "+r"$\bar{n}$", "No PC + Shift + "+r"$n_\text{fid}$"],legend_loc='center right',legend_ncol=1);
    g.export("plot_1d_omegam.pdf")

    g = plots.get_subplot_plotter(width_inch=4)
    print(g,type(g))
    g.plot_1d(chains,"sigma8")
    g.add_x_marker(0.812229,lw=1.5)
    g.settings.legend_fontsize=9
    g.add_legend(["No PC + Shift + "+r"$\bar{n}$", "No PC + Shift + "+r"$n_\text{fid}$"],legend_loc='center left',legend_ncol=1);
    g.export("plot_1d_sigma8.pdf")

#     g = plots.get_subplot_plotter()
#     g.plots_1d([one_pc],["roman_alpha_1"],
#             markers={"roman_alpha_1":0.0},
#             colors=["#E69F00"],
#             ls=['-'],lws=[3])
#     g.export("alpha_1pc.pdf")

#     g = plots.get_subplot_plotter()
#     g.plots_1d([two_pc],["roman_alpha_1","roman_alpha_2"],
#             markers={"roman_alpha_1":0.0,"roman_alpha_2":0.0},
#             colors=["#D55E00"],nx=2,
#             ls=['-'],lws=[3])
#     g.export("alpha_2pc.pdf")

    # g = plots.get_subplot_plotter()
    # g.plots_1d([mcmc11],["roman_alpha_1","roman_alpha_2","roman_alpha_3"],
    #         markers={"roman_alpha_1":0.0,"roman_alpha_2":0.0,"roman_alpha_3":0.0},
    #         colors=["#56B4E9"],nx=3,
    #         ls=['-'],lws=[3])
    # g.export("alpha_3pc.pdf")
    # g.settings.title_limit_fontsize = 14
    # g.triangle_plot([mcmc11
    #                  ],
    #                 ["omegam","sigma8","roman_alpha_1",
    #                  "roman_alpha_2","roman_alpha_3","roman_alpha_4","roman_alpha_5"
    #                 ,"roman_alpha_6","roman_alpha_7","roman_alpha_8","roman_alpha_9","roman_alpha_10"
    #                 ],
    #                  markers={"omegam": 0.3156, "omegab": 0.0492, "sigma8": 0.811, "H0": 67.32,
    #                 "roman_alpha_1": -0.0006333320119154361,
    #                 "roman_alpha_2": 0.002574517895020503,
    #                 "roman_alpha_3": 0.006577603048039325,
    #                 "roman_alpha_4": 0.28971425278177315,
    #                 "roman_alpha_5": -0.004522791743629093,
    #                 "roman_alpha_6": 0.2661174278304146,
    #                 "roman_alpha_7": 0.0015610308292448886,
    #                 "roman_alpha_8": 0.1492990608421329,
    #                 "roman_alpha_9": 0.013365201830143171,
    #                 "roman_alpha_10": -0.2167643230744845,
    #                 }, # planck best-fit 1807.06209
    #                  legend_ncol=1,
    #                  legend_labels=['10 PC'],
    #                  filled=True,
    #                  line_args=[{"ls": "-" ,"color": "C0", "lw":2},]  # no pc + shift)
    # )
    # g.export("triangle_plot_10pc.pdf")
    # return None
roman_pca_mcmc()

# def roman_real_mcmc():
#     mcmc0 = loadMCSamples(path+f"{px}_MCMC0", settings=settings)
#     mcmc1 = loadMCSamples(path+f"{px}_MCMC1", settings=settings)
#     mcmc2 = loadMCSamples(path+f"{px}_MCMC2", settings=settings)

#     g = plots.get_subplot_plotter()
#     g.plots_1d([mcmc0,mcmc1,mcmc2],["omegam","sigma8","w","wa"],
#             markers={"omegam":0.3,"sigma8":0.8120,"w":-1,"wa":0}, # planck best-fit 1807.06209
#             nx=4,legend_ncol=5,
#             legend_labels=['ROMAN_REAL_MCMC0', 'ROMAN_REAL_MCMC1', 'ROMAN_REAL_MCMC2'],
#             colors=['#1b5f6f',"#E69F00","#D55E00"],
#             ls=['-','--','-.'],lws=[3,3,3])
#     g.export("1d_marg.pdf")

#     g = plots.get_subplot_plotter()
#     g.plots_1d([mcmc0,mcmc1,mcmc2],["roman_DZ_S1","roman_DZ_S2","roman_DZ_S3"],
#             markers={"roman_DZ_S1":0.0,"roman_DZ_S2":0.0,"roman_DZ_S3":0.0},
#             colors=['#1b5f6f',"#E69F00","#D55E00"],nx=3,
#             ls=['-','--','-.'],lws=[3,3,3])
#     g.export("1d_marg_alpha_3pc.pdf")
#     return None
# roman_real_mcmc()


# def plot_RMinus1(type_):
#     plt.figure()
#     Rminus1_stop = 0.015
#     max_val = 0
#     progress_0 = f'{px}_MCMC0.progress'
#     progress_1 = f'{px}_MCMC1.progress'
#     progress_2 = f'{px}_MCMC2.progress'
#     progress_3 = f'{px}_MCMC3.progress'

#     # progress = [progress_0,progress_1,progress_2,progress_3]
#     # colors = ['#1b5f6f',"#E69F00","#D55E00","#56B4E9"]
#     # legend_labels=['No PC', '1 PC', '2 PC', '3 PC']
    
#     if type_ == 'roman_pca': 
#         progress = [progress_1,progress_2,progress_3]
#         colors = ["#E69F00","#D55E00","#56B4E9"]
#         legend_labels=['1 PC', '2 PC', '3 PC']

#     elif type_ == 'roman_shift':
#         progress = [progress_0,progress_1,progress_2]
#         colors = ['#1b5f6f',"#E69F00","#D55E00"]
#         legend_labels=['ROMAN_REAL_MCMC0', 'ROMAN_REAL_MCMC1', 'ROMAN_REAL_MCMC2']

#     for p,c,l in zip(progress,colors,legend_labels):
#         Rminus1 = np.genfromtxt(path+p,usecols=(3))
#         Rminus1 = Rminus1[~np.isnan(Rminus1)]
#         vals = list(range(len(Rminus1)))
#         if  np.max(vals) > max_val:
#             max_val = np.max(vals)
#         plt.plot(vals,Rminus1,label=l,marker='o',markersize=3,color=c)

#     plt.hlines(Rminus1_stop,np.min(vals),max_val,color='gray',ls='--',label=f'Rminus_stop: {Rminus1_stop}',lw=3)
#     plt.yscale('log')
#     plt.legend(loc='best',fontsize=13)
#     plt.ylabel('R-1',fontsize=20)
#     plt.xlabel('Points',fontsize=20)
#     plt.xlim(np.min(vals),max_val)
#     plt.ylim(Rminus1_stop-0.1,np.max(Rminus1)+5)
#     plt.savefig('Rminus1.pdf')
#     return None
# # type_ = 'roman_pca'
# type_ = 'roman_shift'
# plot_RMinus1(type_)