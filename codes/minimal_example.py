import numpy as np
import matplotlib.pyplot as plt
from getdist import plots,loadMCSamples,MCSamples

burnin=0.5
settings = {"ignore_rows": burnin}
cocoa_path='/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa'
path = f'{cocoa_path}/cocoa_photoz/results/chains/roman_pca/MCMC'

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

idxes = [57,59,69,79] #[fid,0pc,5pc,10pc]

print(f'USING MCMCs WITH INDEX: {idxes}')

chains=[]
for i in idxes:
    chains.append(loadMCSamples(f"{path}{i}", settings=settings)) 

g = plots.get_subplot_plotter()
g.settings.legend_fontsize = 11

g.triangle_plot(chains,
                ["omegam","sigma8",
                #  markers={"omegam": 0.3, "omegab": 0.04, "sigma8": 0.8277, "H0": 67.32,},
                    ],
                legend_ncol=1,
                markers={"omegam": 0.3, "sigma8": 0.8277},
                legend_labels=["Fiducial","0 PC","5 PCs", "10 PCs"],
                line_args=[{"ls": "-" , "color": "#0072B2", "lw":2.5}, # Fiducial
                        {"ls": ":" , "color": "#E69F00", "lw":2.5},    # 0pc
                        {"ls": "-.", "color": "#009E73", "lw":2.5},    # 5pc
                        {"ls": "--", "color": "#CC79A7", "lw":2.5},    # 10pc
                        ],
                contour_colors=["#0072B2",
                                "#E69F00",
                                "#009E73", 
                                "#CC79A7",
                                ],
                filled=True)
g.export(f"cocoa_photoz/codes/minimal_example/triangle_plot.pdf")

g.triangle_plot(chains[-1],
                [
                    "roman_alpha_1","roman_alpha_2","roman_alpha_3","roman_alpha_4","roman_alpha_5",
                    "roman_alpha_6","roman_alpha_7","roman_alpha_8","roman_alpha_9","roman_alpha_10",
                    ],
                legend_ncol=1,
                legend_labels=["10 PCs"],
                filled=True)
g.export("cocoa_photoz/codes/minimal_example/triangle_plot_10pc.pdf")
