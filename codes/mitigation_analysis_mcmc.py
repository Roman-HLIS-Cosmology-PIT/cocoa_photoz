import numpy as np
import matplotlib.pyplot as plt
from getdist import plots,loadMCSamples,MCSamples

burnin=0.5
settings = {"ignore_rows": burnin}
cocoa_path='/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa'

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

chains_path = f'{cocoa_path}/cocoa_photoz/results/chains/roman_pca'
chains_unw_path = f'{chains_path}/MCMC'
chains_w_path = f'{chains_path}/weighted_pca/MCMC'

"""
# 57: fiducial (zero bias): No PC + No Shift with synthetic xi(nfid) + model xi(fid)
# 59: (highest bias): No PC + No Shift with synthetic xi(nfid) + model xi(nbar)
"""
initial_idx = [57]
initial_chains=[]
for iidx in initial_idx:
    initial_chains.append(loadMCSamples(f"{chains_path}/MCMC{iidx}", settings=settings))

unw_indx  = [61,63,65,67] # 1pc,2pc,3pc,4pc [UN-WEIGHTED PCA]
unw_indx2 = [69,71,73,75] # 5pc,6pc,7pc,8pc [UN-WEIGHTED PCA]
w_indx    = [1,3,5,7]   # 1pc,2pc,3pc,4pc [WEIGHTED PCA]
w_indx2   = [9,11,13,15]   # 5pc,6pc,7pc,8pc [WEIGHTED PCA]

config={"params":["omegam","sigma8"],
        "markers":{"omegam": 0.3, "omegab": 0.04, "sigma8": 0.8277, "H0": 67.32,},
        "legend_ncol":1,
        "line_args":[{"ls": "-", "color": "gray", "lw":1.5},
                     {"ls": "-" , "color": "#0072B2", "lw":2.5},    
                     {"ls": ":" , "color": "#E69F00", "lw":2.5},  
                     {"ls": "--", "color": "#009E73", "lw":2.5},  
                     {"ls": "-.",  "color": "#CC79A7", "lw":2.5 }],
        "contour_colors":["gray","#0072B2","#E69F00","#009E73", "#CC79A7",],
        "filled":True}


#############################################
#############################################
def func1():
    chains_unw_pca = initial_chains.copy()

    label_unw = "57,"

    for i in unw_indx:
        chains_unw_pca.append(loadMCSamples(f"{chains_unw_path}{i}", settings=settings)) 
        label_unw+=f"{i},"

    g = plots.get_subplot_plotter()
    g.settings.legend_fontsize = 11
    g.triangle_plot(roots=chains_unw_pca,**config,legend_labels=["Fiducial","1 PC","2 PC", "3 PC", "4 PC"])
    g.export(f"triangle_plot_unweighted_{label_unw}.pdf")    
    return None
#############################################
#############################################
def func2():
    chains_unw_pca = initial_chains.copy()

    label_unw = "57,"

    for i in unw_indx2:
        chains_unw_pca.append(loadMCSamples(f"{chains_unw_path}{i}", settings=settings)) 
        label_unw+=f"{i},"

    g = plots.get_subplot_plotter()
    g.settings.legend_fontsize = 11
    g.triangle_plot(roots=chains_unw_pca,**config,legend_labels=["Fiducial","5 PC","6 PC", "7 PC", "8 PC"])
    g.export(f"triangle_plot_unweighted_{label_unw}.pdf")
    return None
#############################################
#############################################
def func3():
    chains_w_pca = initial_chains.copy()

    label_w = "57,"

    for i in w_indx:
        chains_w_pca.append(loadMCSamples(f"{chains_w_path}{i}", settings=settings)) 
        label_w+=f"{i},"

    g = plots.get_subplot_plotter()
    g.settings.legend_fontsize = 11
    g.triangle_plot(roots=chains_w_pca,**config,legend_labels=["Fiducial","1 PC","2 PC", "3 PC", "4 PC"])
    g.export(f"triangle_plot_weighted_{label_w}.pdf")
    return None
#############################################
#############################################
def func4():
    chains_w_pca = initial_chains.copy()

    label_w = "57,"

    for i in w_indx2:
        chains_w_pca.append(loadMCSamples(f"{chains_w_path}{i}", settings=settings)) 
        label_w+=f"{i},"

    g = plots.get_subplot_plotter()
    g.settings.legend_fontsize = 11
    g.triangle_plot(roots=chains_w_pca,**config,legend_labels=["Fiducial","5 PC","6 PC", "7 PC", "8 PC"])
    g.export(f"triangle_plot_weighted_{label_w}.pdf")
    return None

# func4()
