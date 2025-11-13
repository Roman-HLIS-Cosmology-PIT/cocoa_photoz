import numpy as np
import matplotlib.pyplot as plt
from getdist import plots,loadMCSamples,MCSamples
import re
import inspect
import argparse


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

colors1 = [
    "#1b5f6f",  # dark teal
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#882255",  # wine
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#999999",  # grey
    "#117733",  # dark green
]


########################################################
import matplotlib
plt.rcParams['figure.figsize'] = (3.5, 2.5)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['xtick.bottom'] = True
matplotlib.rcParams['xtick.top'] = False
matplotlib.rcParams['ytick.right'] = False
matplotlib.rcParams['axes.edgecolor'] = 'black'
matplotlib.rcParams['axes.linewidth'] = '1.0'
matplotlib.rcParams['axes.labelsize'] = 'medium'
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['grid.linewidth'] = '0.0'
matplotlib.rcParams['grid.alpha'] = '0.18'
matplotlib.rcParams['grid.color'] = 'lightgray'
matplotlib.rcParams['legend.labelspacing'] = 0.77
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.format'] = 'pdf'
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['font.size'] = 15
########################################################

chains_path = f'{cocoa_path}/cocoa_photoz/results/chains/roman_pca'
chains_unw_path = f'{chains_path}/MCMC'
chains_w_path = f'{chains_path}/weighted_pca/MCMC'
chains_sc1bd5_x_sc1bd4_path = f'{chains_path}/weighted_pca_Model_sc1bd5_Fiducial_sc1bd4/MCMC'
chains_sc1bd6_x_sc1bd4_path = f'{chains_path}/weighted_pca_Model_sc1bd6_Fiducial_sc1bd4/MCMC'
chains_sc1bd7_x_sc1bd4_path = f'{chains_path}/weighted_pca_Model_sc1bd7_Fiducial_sc1bd4/MCMC'
chain_baseline = f'{chains_path}/weighted_pca/MCMC_BASELINE'

"""
# 57: fiducial (zero bias): No PC + No Shift with synthetic xi(nfid) + model xi(fid)
# 59: (highest bias): No PC + No Shift with synthetic xi(nfid) + model xi(nbar)
"""
initial_idx = [57]
initial_chains=[]
for iidx in initial_idx:
    samples = loadMCSamples(f"{chains_path}/MCMC{iidx}", settings=settings)
    sigma8 = samples["sigma8"]
    omegam = samples["omegam"]
    S8 = sigma8*np.sqrt(omegam/0.3)
    samples.addDerived(S8,name="S8",label="S_8")
    initial_chains.append(samples)

### NON MIXING SCENARIOS ###
unw_indx  = [61,63,65,67]        ##  1pc,2pc,3pc,4pc    [UN-WEIGHTED PCA]

w_indx_le = [81,82,84,86]   ## 1pc,2pc,3pc,4pc,5pc  [WEIGHTED PCA] (least extreme - r = 3233)
w_indx_in = [102,103,105,107] ## 1pc,2pc,3pc,4pc,5pc  [WEIGHTED PCA] (intermediate - r = 0)
w_indx_me = [123,124,126,128] ## 1pc,2pc,3pc,4pc,5pc  [WEIGHTED PCA] (most extreme - r = 869)

w_indx_le = [81,84,86]#list(range(81,87))   ## 0pc,1pc,2pc,3pc,4pc,5pc  [WEIGHTED PCA] (least extreme - r = 3233)
w_indx_in = [102,105,107]#list(range(102,108)) ## 0pc,1pc,2pc,3pc,4pc,5pc  [WEIGHTED PCA] (intermediate - r = 0)
w_indx_me = [123,126,128]#list(range(123,129)) ## 0pc,1pc,2pc,3pc,4pc,5pc  [WEIGHTED PCA] (most extreme - r = 869)

### MIXING SCENARIOS ###
sc1bd5_x_sc1bd4_indx = list(range(6))
sc1bd6_x_sc1bd4_indx = list(range(6))
sc1bd7_x_sc1bd4_indx = list(range(6))

combined_indx = [81,102,123]
baseline  = [22,23,24]

bconfig = {"params":["omegam","S8"],
        "markers":{"omegam": 0.3, "S8":0.8277*(0.3/0.3)**0.5},
        "legend_ncol":1,
        "line_args":[{"ls": "-" , "color": colors1[0], "lw":1.5},
                     {"ls": "--", "color": colors1[1], "lw":1.5},
                     {"ls": "-.", "color": colors1[2], "lw":1.5},
                     {"ls": ":" , "color": colors1[3], "lw":1.5},
                     {"ls": "-" , "color": colors1[4], "lw":1.5},
                     {"ls": "--", "color": colors1[5], "lw":1.5},
                     ],
        "contour_colors":[colors1[0],colors1[1],colors1[2],colors1[3],colors1[4],colors1[5]],
        "contour_ls":["-","--","-.",":","-","--"],
        "contour_lws":[1.5,1.5,1.5,1.5,1.5,1.5],
        "filled":False}

config = lambda c: {"params":["omegam","S8"],
        "markers":{"omegam": 0.3, "S8":0.8277*(0.3/0.3)**0.5},
        "legend_ncol":1,
        "line_args":[{"ls": "-" , "color": c, "lw":1.5},
                     {"ls": "--", "color": c, "lw":1.5},
                     {"ls": ":" , "color": c, "lw":1.5},
                     ],
        "contour_colors":[c,c,c],
        "contour_ls":["-","--",":"],
        "contour_lws":[1.5,1.5,1.5],
        "legend_labels":["0 PC","3 PC","5 PC"],
        "filled":False}

#############################################
#############################################

def fig_to_array(fig, dpi=600):
    from io import BytesIO
    """Convert a Matplotlib Figure to a numpy RGBA array."""
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight', transparent=True)
    buf.seek(0)
    img = plt.imread(buf)
    buf.close()
    return img

def triangle_plot():
    chains_le = [] #initial_chains.copy()
    label_le=''
    for i in w_indx_le:
        samples = loadMCSamples(f"{chains_w_path}{i}", settings=settings)
        sigma8,omegam = samples["sigma8"],samples["omegam"]
        S8 = sigma8*np.sqrt(omegam/0.3)
        samples.addDerived(S8,name="S8",label="S_8")
        chains_le.append(samples) 
        label_le+=f"{i},"
    g_le = plots.get_subplot_plotter()
    g_le.settings.legend_fontsize = 11
    g_le.triangle_plot(roots=chains_le,**config("#1b5f6f"))
    fig_name=f"poster_triangle_plot_weighted_{label_le}_Model_sc1bd4_Fiducial_sc1bd4.pdf"
    g_le.export(fig_name)
    print(fig_name)

    chains_in = [] #initial_chains.copy()
    label_in=''
    for i in w_indx_in:
        samples = loadMCSamples(f"{chains_w_path}{i}", settings=settings)
        sigma8,omegam = samples["sigma8"],samples["omegam"]
        S8 = sigma8*np.sqrt(omegam/0.3)
        samples.addDerived(S8,name="S8",label="S_8")
        chains_in.append(samples) 
        label_in+=f"{i},"
    g_in = plots.get_subplot_plotter()
    g_in.settings.legend_fontsize = 11
    g_in.triangle_plot(roots=chains_in,**config("#E69F00"))
    fig_name=f"poster_triangle_plot_weighted_{label_in}_Model_sc1bd4_Fiducial_sc1bd4.pdf"
    g_in.export(fig_name)
    print(fig_name)

    chains_me = [] #initial_chains.copy()
    label_me=''
    for i in w_indx_me:
        samples = loadMCSamples(f"{chains_w_path}{i}", settings=settings)
        sigma8,omegam = samples["sigma8"],samples["omegam"]
        S8 = sigma8*np.sqrt(omegam/0.3)
        samples.addDerived(S8,name="S8",label="S_8")
        chains_me.append(samples) 
        label_me+=f"{i},"
    g_me = plots.get_subplot_plotter()
    g_me.settings.legend_fontsize = 11
    g_me.triangle_plot(roots=chains_me,**config("#882255"))
    fig_name=f"poster_triangle_plot_weighted_{label_me}_Model_sc1bd4_Fiducial_sc1bd4.pdf"
    g_me.export(fig_name)
    print(fig_name)
    return None

def plot_chisq():
    colors1 = [
        "#1b5f6f",  # dark teal
        "#E69F00",  # orange
        "#56B4E9",  # sky blue
        "#882255",  # wine
        "#009E73",  # bluish green
        "#D55E00",  # vermillion
        "#F0E442",  # yellow
        "#0072B2",  # blue
        "#CC79A7",  # reddish purple
        "#999999",  # grey
        "#117733",  # dark green
    ]
    indexes = [123,124,125,126,127,128]
    config = {"param":"chi2",
              "legend_ncol":1,
              "line_args":[
                  {"ls": "-" , "color": colors1[0], "lw":1.5},
                  {"ls": "--" , "color": colors1[1], "lw":1.5},
                  {"ls": "-." , "color": colors1[2], "lw":1.5},
                  {"ls": ":" , "color": colors1[3], "lw":1.5},
                  {"ls": "-" , "color": colors1[4], "lw":1.5},
                  {"ls": "--" , "color": colors1[5], "lw":1.5},
                  {"ls": "--" , "color": "black", "lw":2.5},
                  ],
            #   "contour_lws":[1.5,1.5,1.5,1.5,1.5,1.5,1.5],
              "legend_labels":["0 PC","1 PC","2 PC","3 PC","4 PC","5 PC","Shift"],
            }
    chains_me = []
    label_me=''
    for i in indexes:
        print(f"chain index {i}")
        samples = loadMCSamples(f"{chains_w_path}{i}", settings=settings)
        chains_me.append(samples) 
        label_me+=f"{i},"
    for i in [24]:
        print(f"chain index {i}")
        samples = loadMCSamples(f"{chain_baseline}{i}", settings=settings)
        chains_me.append(samples) 
        label_me+=f"{i},"
    plt.figure(figsize=(10, 10))    
    g_me = plots.get_single_plotter(width_inch=6)
    g_me.settings.legend_fontsize = 11
    g_me.plot_1d(roots=chains_me,**config)
    plt.legend(config["legend_labels"],fontsize=10)
    plt.title("Most Extreme")
    fig_name=f"plot1d_chisq_{label_me}_Model_sc1bd4_Fiducial_sc1bd4.pdf"
    g_me.export(fig_name)
    print(fig_name)
    return None

def plot_2d():
    chains = [] 
    label = ''
    for i in combined_indx:
        samples = loadMCSamples(f"{chains_w_path}{i}", settings=settings)
        sigma8,omegam = samples["sigma8"],samples["omegam"]
        S8 = sigma8*np.sqrt(omegam/0.3)
        samples.addDerived(S8,name="S8",label="S_8")
        chains.append(samples) 
        label+=f"{i},"
    g = plots.get_single_plotter()
    g.settings.legend_fontsize = 11
    g.plot_2d(chains, "omegam", "S8")
    
    fig_name=f"plot_2d_{label}.pdf"
    g.export(fig_name)
    print(fig_name)
    return None

def get_inline_latex_lcdm():
    print(f"EXECUTING get_inline_latex_lcdm @ LINE {inspect.currentframe().f_lineno}")
    NPCs = 5
    combined_indx = list(range(81,81+NPCs+1)) + \
                    list(range(102,102+NPCs+1)) + \
                    list(range(123,123+NPCs+1))
    baseline  = [22,23,24] 
    baseline_extended  = [22,23,24,  # Prior: 0.003
                          25,26,27,  # Prior: 0.01
                        #   28,29,30 # Prior: 0.1
                          ]  

    chains = [] 
    chains_baseline = [] 
    label = ''
    pattern = r"=\s*([\d.]+)\^\{\+([\d.]+)\}_\{\-([\d.]+)\}"

    ###################
    ### LOAD CHAINS ###
    ###################
    for i in combined_indx:
        samples = loadMCSamples(f"{chains_w_path}{i}", settings=settings)
        sigma8,omegam = samples["sigma8"],samples["omegam"]
        S8 = sigma8*np.sqrt(omegam/0.3)
        samples.addDerived(S8,name="S8",label="S_8")
        chains.append(samples) 

    for i in baseline_extended:
        samples = loadMCSamples(f"{chain_baseline}{i}", settings=settings)
        sigma8,omegam = samples["sigma8"],samples["omegam"]
        S8 = sigma8*np.sqrt(omegam/0.3)
        samples.addDerived(S8,name="S8",label="S_8")
        chains_baseline.append(samples) 

    ######################
    ### PCA 95% ERRORS ###
    ######################
    omegams_pca = []
    S8s_pca = []
    for ch in chains:
        omegam_stats = ch.getInlineLatex("omegam", limit=2)
        S8_stats = ch.getInlineLatex("S8", limit=2)
        match1 = re.search(pattern, omegam_stats)
        match2 = re.search(pattern, S8_stats)
        print("PCs:",omegam_stats,S8_stats)  
        if match1:
            central = float(match1.group(1))
            upper = float(match1.group(2))
            lower = -float(match1.group(3))  # make it negative explicitly
            omegams_pca.append([central,lower,upper])
        if match2:
            central = float(match2.group(1))
            upper = float(match2.group(2))
            lower = -float(match2.group(3))  # make it negative explicitly
            S8s_pca.append([central,lower,upper])

    ###################################
    ### BASELINE (SHIFT) 95% ERRORS ###
    ###################################
    omegams_baseline = []
    S8s_baseline = []
    for ch in chains_baseline:
        omegam_stats = ch.getInlineLatex("omegam", limit=2)
        S8_stats = ch.getInlineLatex("S8", limit=2)
        match1 = re.search(pattern, omegam_stats)
        match2 = re.search(pattern, S8_stats)
        print("Baseline",omegam_stats,S8_stats)                    
        if match1:
            central = float(match1.group(1))
            upper = float(match1.group(2))
            lower = -float(match1.group(3))  # make it negative explicitly
            omegams_baseline.append([central,lower,upper])
        if match2:
            central = float(match2.group(1))
            upper = float(match2.group(2))
            lower = -float(match2.group(3))  # make it negative explicitly
            S8s_baseline.append([central,lower,upper])
    
    print("omegams_pca:",omegams_pca)
    print("S8s_pca:",S8s_pca)
    print("omegams_baseline:",omegams_baseline)
    print("S8s_baseline:",S8s_baseline)
    
    omegam_fid = 0.3
    S8_fid = 0.8277*(0.3/0.3)**0.5

    ############
    ### PLOT ###
    ############
    fig = plt.figure(figsize=(10, 4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, wspace=0.0, hspace=0.0, width_ratios=[3,1])

    # Left column axes (share y)
    axes00 = fig.add_subplot(gs[0,0])
    axes10 = fig.add_subplot(gs[1,0],sharex=axes00)  # shares y with top-left

    # Right column axes (share y)
    axes01 = fig.add_subplot(gs[0,1],sharey=axes00)
    axes11 = fig.add_subplot(gs[1,1],sharey=axes10)  # shares y with top-right

    axes00.axhline(y=omegam_fid, color='gray', linestyle='--', linewidth=1.2)
    axes10.axhline(y=S8_fid, color='gray', linestyle='--', linewidth=1.2)
    axes01.axhline(y=omegam_fid, color='gray', linestyle='--', linewidth=1.2)
    axes11.axhline(y=S8_fid, color='gray', linestyle='--', linewidth=1.2)

    ecolors = ["#1b5f6f"]*(NPCs+1)+\
              ["#E69F00"]*(NPCs+1)+\
              ["#882255"]*(NPCs+1)
    
    fmts = ["o"]*(NPCs+1)+\
           ["^"]*(NPCs+1)+\
           ["s"]*(NPCs+1)
    
    dx = [-.1]*(NPCs+1)+\
         [0]*(NPCs+1)+\
         [.1]*(NPCs+1)

    cases = list(range(len(combined_indx)))
    npcs = list(range(NPCs+1))
    pairs = [(npcs[i % len(npcs)], c) for i, c in enumerate(cases)]
    
    for p in pairs:
        npc = p[0]
        case = p[1]
        axes00.errorbar(x=npc+dx[case], y=omegams_pca[case][0], yerr=[[abs(omegams_pca[case][1])], [omegams_pca[case][2]]],
                        fmt=fmts[case], capsize=4, color=ecolors[case],ecolor=ecolors[case], elinewidth=1.2)
        axes10.errorbar(x=npc+dx[case], y=S8s_pca[case][0], yerr=[[abs(S8s_pca[case][1])], [S8s_pca[case][2]]],
                        fmt=fmts[case], capsize=4, color=ecolors[case], ecolor=ecolors[case], elinewidth=1.2)
    
    ecolors = ["#1b5f6f","#E69F00","#882255", # Prior 0.003  
               "#1b5f6f","#E69F00","#882255", # Prior 0.01  
               ]
    mfcs = ["#1b5f6f","#E69F00","#882255", # Prior 0.003  
             "none","none","none",         # Prior 0.01  
               ]
    fmts = ["o","^","s", # Prior 0.003  
            "o","^","s", # Prior 0.01  
            ]
    dx = [-0.01,0.00,0.01,    # Prior 0.003  
          -0.008,0.002,0.012  # Prior 0.01  
          ]
    for case in [0,1,2, 
                 3,4,5, 
                 ]:
        axes01.errorbar(x=dx[case], y=omegams_baseline[case][0], yerr=[[abs(omegams_baseline[case][1])], [omegams_baseline[case][2]]],
                    fmt=fmts[case], capsize=4, color=ecolors[case],ecolor=ecolors[case],mfc=mfcs[case], elinewidth=1.2)
        axes11.errorbar(x=dx[case], y=S8s_baseline[case][0], yerr=[[abs(S8s_baseline[case][1])], [S8s_baseline[case][2]]],
                    fmt=fmts[case], capsize=4, color=ecolors[case],ecolor=ecolors[case],mfc=mfcs[case], elinewidth=1.2)
    
    print(omegams_pca)
    print(S8s_pca)
    
    axes00.set_ylabel(r"$\Omega_m$")
    axes10.set_ylabel(r"$S_8$")
    axes10.set_xlabel(r"# of PCs")
    axes01.set_title(r"Shift $\Delta_z^i$")
    axes10.set_xticks([0, 1, 2, 3, 4, 5])
    axes10.set_xlim(-0.2,5.2)
    axes11.set_xlabel(r"Prior")
    axes11.set_xticks(list(np.sort(dx)))
    axes11.set_xticklabels(['a', 'b']*3)
    axes00.xaxis.set_visible(False)
    # for ax in [axes01, axes11]:
    #     ax.yaxis.set_visible(False)
    #     ax.xaxis.set_visible(False)    
    for ax in [axes01]:
        ax.yaxis.set_visible(False)
        ax.xaxis.set_visible(False)
    for ax in [axes11]:
        ax.yaxis.set_visible(False)    

    # Legend for marker types
    from matplotlib.lines import Line2D
    marker_labels = {'o': 'Least Extreme', '^': 'Intermediate', 's': 'Most Extreme'}
    colors = ["#1b5f6f", "#E69F00", "#882255"]

    legend_handles = [
    Line2D([0], [0], marker=m, color=c, linestyle='None', markersize=8, label=l)
    for (m, l), c in zip(marker_labels.items(), colors)]

    # Place legend above axes[0]
    axes00.legend(
    handles=legend_handles,
    ncol=3,               # horizontal layout
    title=None,           # no title
    fontsize=13,
    loc='lower center',   # anchor at bottom center of bbox_to_anchor
    bbox_to_anchor=(0.5, 1.0)  # x=0.5 centered, y=1.05 slightly above axes
    )

    fig_name="mitigation_least_intermediate_most_extreme.pdf"
    plt.savefig(fig_name) # Paper Figure 11 
    print(fig_name)
    return None

def bias_error_mitigation_lcdm(submit_bias=False,want_plot=True,submit_errors=False):
    """
    submit_bias: submit job in case you don't have the file samps_stats.txt
    want_plot: if you have the file, then you likely want to make the plot, otherwise, not. 
    (Likely need to run twice if the file doesn't exist.)
    submit_errors: if true, plot 0.5sigma errors (see arxiv 1809.01146 Figure 11) 
    """
    print(f"EXECUTING bias_error_mitigation_lcdm @ LINE {inspect.currentframe().f_lineno}")
    NPCs = 5
    combined_indx = list(range(81,81+NPCs+1)) + \
                    list(range(102,102+NPCs+1)) + \
                    list(range(123,123+NPCs+1))

                         # least, intermediate, most - extreme
    baseline_extended  = [22,23,24,  # Prior: 0.003
                          25,26,27,  # Prior: 0.01
                          ]  

    #######################################################
    #### IF THE FILE DON'T EXIST, UNCOMMENT AND SUBMIT ####
    #######################################################
    if submit_bias:
        print("PCA STARTS")
        for i in combined_indx:
            samples = loadMCSamples(f"{chains_w_path}{i}", settings=settings)
            samps_stats = samples.getLikeStats()

            print(f"sample {i}")
            print(samps_stats)
        print("PCA ENDS")

        print("BASELINE STARTS")
        for i in baseline_extended:
            samples_base = loadMCSamples(f"{chain_baseline}{i}", settings=settings)
            samps_stats_base = samples_base.getLikeStats()

            print(f"sample {i}")
            print(samps_stats_base)
        print("BASELINE ENDS")
        print('DONE ALL - PCA AND BASELINE')    
    #######################################################
    #### IF THE FILE DON'T EXIST, UNCOMMENT AND SUBMIT #### 
    #######################################################  
    sigma_stas=2
    if submit_errors:
        print("PCA AND BASELINE 0.5*SIGMA ERRORS - STARTS")
        """include_errors: if true, plot 0.5sigma errors (see arxiv 1809.01146 Figure 11) """
        chains=[]
        chains_baseline=[]
        pattern_asymmetric = r"=\s*([\d.]+)\^\{\+([\d.]+)\}_\{\-([\d.]+)\}"
        pattern_symmetric = r"(?:\\pm|±)\s*([\d.]+)"
        ###################
        ### LOAD CHAINS ###
        ###################
        for i in combined_indx:
            samples = loadMCSamples(f"{chains_w_path}{i}", settings=settings)
            sigma8,omegam = samples["sigma8"],samples["omegam"]
            S8 = sigma8*np.sqrt(omegam/0.3)
            samples.addDerived(S8,name="S8",label="S_8")
            chains.append(samples) 

        for i in baseline_extended:
            samples = loadMCSamples(f"{chain_baseline}{i}", settings=settings)
            sigma8,omegam = samples["sigma8"],samples["omegam"]
            S8 = sigma8*np.sqrt(omegam/0.3)
            samples.addDerived(S8,name="S8",label="S_8")
            chains_baseline.append(samples) 

        ######################
        ### PCA 95% ERRORS ###
        ######################
        omegams_pca_errors = []
        S8s_pca_errors = []
        for ch in chains:
            omegam_stats = ch.getInlineLatex("omegam", limit=sigma_stas)
            S8_stats = ch.getInlineLatex("S8", limit=sigma_stas)
            match1_asy = re.search(pattern_asymmetric, omegam_stats)
            match2_asy = re.search(pattern_asymmetric, S8_stats)
            match1_sy = re.search(pattern_symmetric, omegam_stats)
            match2_sy = re.search(pattern_symmetric, S8_stats)
            print("PCs:",omegam_stats,S8_stats)    
            # ASSYMETRIC OR SYMMETRIC ERRORS - OMEGA       
            if match1_asy is not None:
                upper, lower = float(match1_asy.group(2)), float(match1_asy.group(3))  
                half_sigma = (upper+lower)*0.5
                omegams_pca_errors.append(half_sigma)
                print('0.5sigma omega [asy]: ', half_sigma)
            if match1_sy is not None:
                print("match1_sy:",match1_sy)
                upper = float(match1_sy.group(1))
                omegams_pca_errors.append(upper)
                print('0.5sigma omega [sy]: ', upper)
            # ASSYMETRIC OR SYMMETRIC ERRORS - S8       
            if match2_asy is not None:
                upper, lower = float(match2_asy.group(2)), float(match2_asy.group(3)) 
                half_sigma = (upper+lower)*0.5
                S8s_pca_errors.append(half_sigma)
                print('0.5sigma S8 [asy]: ', half_sigma)
            if match2_sy is not None:
                print("match2_sy:",match2_sy)
                upper = float(match2_sy.group(1))
                S8s_pca_errors.append(upper)
                print('0.5sigma S8 [sy]: ', upper)
        print("pca omegam errors (0.5sigma): ", omegams_pca_errors)        
        print("pca S8 errors (0.5sigma): ", S8s_pca_errors)        
        print("PCA ENDS")
        ###################################
        ### BASELINE (SHIFT) 95% ERRORS ###
        ###################################
        omegams_baseline_errors = []
        S8s_baseline_errors = []
        for ch in chains_baseline:
            omegam_stats = ch.getInlineLatex("omegam", limit=sigma_stas)
            S8_stats = ch.getInlineLatex("S8", limit=sigma_stas)
            match1_asy = re.search(pattern_asymmetric, omegam_stats)
            match2_asy = re.search(pattern_asymmetric, S8_stats)
            match1_sy = re.search(pattern_symmetric, omegam_stats)
            match2_sy = re.search(pattern_symmetric, S8_stats)
            print("Baseline",omegam_stats,S8_stats)  
            # ASSYMETRIC OR SYMMETRIC ERRORS - OMEGA                  
            if match1_asy is not None:
                upper, lower = float(match1_asy.group(2)), float(match1_asy.group(3)) 
                half_sigma = (upper+lower)*0.5
                omegams_baseline_errors.append(half_sigma)
            if match1_sy is not None:
                upper = float(match1_sy.group(1))
                omegams_baseline_errors.append(upper)
            # ASSYMETRIC OR SYMMETRIC ERRORS - S8                  
            if match2_asy is not None:
                upper, lower = float(match2_asy.group(2)), float(match2_asy.group(3))
                half_sigma = (upper+lower)*0.5
                S8s_baseline_errors.append(half_sigma)
            if match2_sy is not None:
                upper = float(match2_sy.group(1))
                S8s_baseline_errors.append(upper)

        print("baseline omegam errors (0.5sigma): ", omegams_baseline_errors)        
        print("baseline S8 errors (0.5sigma): ", S8s_baseline_errors)  
        print("BASELINE ENDS")
        print("PCA AND BASELINE 0.5*SIGMA ERRORS - ENDS")
    if want_plot:

        # Fiducial values
        omegam_fid = 0.3
        S8_fid = 0.8277*(0.3/0.3)**0.5

        # Empty list for bias
        omegams_pca_bias = []
        S8s_pca_bias = []

        ############################
        ### LOAD BEST FIT VALUES ###
        ############################
        with open("/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa/samps_stats_bestfit.txt","r") as f:
            lines = f.readlines()  
        for line in lines:
            if re.match(r'^omegam\b', line):
                parts = line.split()
                omegam_bestfit = float(parts[1])
                omegams_pca_bias.append(abs(omegam_fid-omegam_bestfit))  
            if re.match(r'^sigma8*\b', line):
                # JUST BECAUSE omegam COME FIRST!
                parts = line.split()
                sigma8_bestfit = float(parts[1])
                S8_bestfit = sigma8_bestfit*np.sqrt(omegam_bestfit/0.3)
                S8s_pca_bias.append(abs(S8_fid-S8_bestfit))  

        #################################
        ### LOAD BEST 0.5SIGMA ERRORS ###
        #################################
        import ast
        print(f"USING {sigma_stas} SIGMAS")
        with open(f"/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa/samps_stats_0.5_{sigma_stas}sigma_erros.txt","r") as f:
            lines = f.readlines()  
        for line in lines:
            # PCA
            if line.startswith("pca omegam errors (0.5sigma):"):
                _, list_str = line.split(":", 1)
                omegams_pca_errors = ast.literal_eval(list_str.strip())
            if line.startswith("pca S8 errors (0.5sigma):"):
                _, list_str = line.split(":", 1)
                S8s_pca_errors = ast.literal_eval(list_str.strip())
            # BASELINE
            if line.startswith("baseline omegam errors (0.5sigma):"):
                _, list_str = line.split(":", 1)
                omegams_baseline_errors = ast.literal_eval(list_str.strip())
            if line.startswith("baseline S8 errors (0.5sigma):"):
                _, list_str = line.split(":", 1)
                S8s_baseline_errors = ast.literal_eval(list_str.strip())

        ############
        ### PLOT ###
        ############
        fig = plt.figure(figsize=(10, 4), constrained_layout=True)
        gs = fig.add_gridspec(2, 2, wspace=0.0, hspace=0.0, width_ratios=[3,1])

        # Left column axes (share y)
        axes00 = fig.add_subplot(gs[0,0])
        axes10 = fig.add_subplot(gs[1,0],sharex=axes00)  # shares y with top-left

        # Right column axes (share y)
        axes01 = fig.add_subplot(gs[0,1],sharey=axes00)
        axes11 = fig.add_subplot(gs[1,1],sharey=axes10,sharex=axes01)  # shares y with top-right

        colors = ["#1b5f6f","#E69F00","#882255"]
        fmts = ["o","^","s"]
        labels=['Least Extreme', 'Intermediate', 'Most Extreme']

        ################
        ### PCA BIAS ###
        ################
        for k in [0,1,2]:
            axes00.plot(list(range(len(omegams_pca_bias[6*k:6*(k+1)]))),omegams_pca_bias[6*k:6*(k+1)],
                    color=colors[k], marker=fmts[k], label=labels[k])
        for k in [0,1,2]:
            axes10.plot(list(range(len(S8s_pca_bias[6*k:6*(k+1)]))),S8s_pca_bias[6*k:6*(k+1)],
                    color=colors[k], marker=fmts[k], label=labels[k])
        ##################
        ### PCA ERRORS ###
        ##################
        print("omegams_pca_errors: ",omegams_pca_errors)
        print("S8s_pca_errors: ",S8s_pca_errors)
        for k in [0,1,2]:
            axes00.plot(list(range(len(omegams_pca_errors[6*k:6*(k+1)]))),omegams_pca_errors[6*k:6*(k+1)],
                    color=colors[k], marker=fmts[k], label=labels[k],ls="--",alpha=0.5)
        for k in [0,1,2]:
            axes10.plot(list(range(len(S8s_pca_errors[6*k:6*(k+1)]))),S8s_pca_errors[6*k:6*(k+1)],
                    color=colors[k], marker=fmts[k], label=labels[k],ls="--",alpha=0.5)
        
        #############################
        ### BASELINE (SHIFT) BIAS ###
        #############################
        # Baseline values is on the same list as pca.
        # Lines
        axes01.plot([0,1],[omegams_pca_bias[18],omegams_pca_bias[21]],color=colors[0]) # Least extreme 0.003 & 0.01
        axes01.plot([0,1],[omegams_pca_bias[19],omegams_pca_bias[22]],color=colors[1]) # Intermediate  0.003 & 0.01
        axes01.plot([0,1],[omegams_pca_bias[20],omegams_pca_bias[23]],color=colors[2]) # Most extreme  0.003 & 0.01
        axes11.plot([0,1],[S8s_pca_bias[18],S8s_pca_bias[21]]        ,color=colors[0]) # Least extreme 0.003 & 0.01
        axes11.plot([0,1],[S8s_pca_bias[19],S8s_pca_bias[22]]        ,color=colors[1]) # Intermediate  0.003 & 0.01
        axes11.plot([0,1],[S8s_pca_bias[20],S8s_pca_bias[23]]        ,color=colors[2]) # Most extreme  0.003 & 0.01
        # Dots
        axes01.scatter([0,1],[omegams_pca_bias[18],omegams_pca_bias[21]],color=[colors[0],colors[0]],marker=fmts[0],facecolors=[colors[0],'none']) # Least extreme 0.003 & 0.01
        axes01.scatter([0,1],[omegams_pca_bias[19],omegams_pca_bias[22]],color=[colors[1],colors[1]],marker=fmts[1],facecolors=[colors[1],'none']) # Intermediate  0.003 & 0.01
        axes01.scatter([0,1],[omegams_pca_bias[20],omegams_pca_bias[23]],color=[colors[2],colors[2]],marker=fmts[2],facecolors=[colors[2],'none']) # Most extreme  0.003 & 0.01
        axes11.scatter([0,1],[S8s_pca_bias[18],S8s_pca_bias[21]]        ,color=[colors[0],colors[0]],marker=fmts[0],facecolors=[colors[0],'none']) # Least extreme 0.003 & 0.01
        axes11.scatter([0,1],[S8s_pca_bias[19],S8s_pca_bias[22]]        ,color=[colors[1],colors[1]],marker=fmts[1],facecolors=[colors[1],'none']) # Intermediate  0.003 & 0.01
        axes11.scatter([0,1],[S8s_pca_bias[20],S8s_pca_bias[23]]        ,color=[colors[2],colors[2]],marker=fmts[2],facecolors=[colors[2],'none']) # Most extreme  0.003 & 0.01

        ###############################
        ### BASELINE (SHIFT) ERRORS ###
        ###############################
        # Baseline values is on the same list as pca.
        # Lines
        axes01.plot([0,1],[omegams_pca_errors[0],omegams_pca_errors[3]],color=colors[0],ls="--",alpha=0.5) # Least extreme 0.003 & 0.01
        axes01.plot([0,1],[omegams_pca_errors[1],omegams_pca_errors[4]],color=colors[1],ls="--",alpha=0.5) # Intermediate  0.003 & 0.01
        axes01.plot([0,1],[omegams_pca_errors[2],omegams_pca_errors[5]],color=colors[2],ls="--",alpha=0.5) # Most extreme  0.003 & 0.01
        axes11.plot([0,1],[S8s_pca_errors[0],S8s_pca_errors[3]]        ,color=colors[0],ls="--",alpha=0.5) # Least extreme 0.003 & 0.01
        axes11.plot([0,1],[S8s_pca_errors[1],S8s_pca_errors[4]]        ,color=colors[1],ls="--",alpha=0.5) # Intermediate  0.003 & 0.01
        axes11.plot([0,1],[S8s_pca_errors[2],S8s_pca_errors[5]]        ,color=colors[2],ls="--",alpha=0.5) # Most extreme  0.003 & 0.01
        # Dots
        axes01.scatter([0,1],[omegams_pca_errors[0],omegams_pca_errors[3]],color=[colors[0],colors[0]],marker=fmts[0],facecolors=[colors[0],'none'],alpha=0.5) # Least extreme 0.003 & 0.01
        axes01.scatter([0,1],[omegams_pca_errors[1],omegams_pca_errors[4]],color=[colors[1],colors[1]],marker=fmts[1],facecolors=[colors[1],'none'],alpha=0.5) # Intermediate  0.003 & 0.01
        axes01.scatter([0,1],[omegams_pca_errors[2],omegams_pca_errors[5]],color=[colors[2],colors[2]],marker=fmts[2],facecolors=[colors[2],'none'],alpha=0.5) # Most extreme  0.003 & 0.01
        axes11.scatter([0,1],[S8s_pca_errors[0],S8s_pca_errors[3]]        ,color=[colors[0],colors[0]],marker=fmts[0],facecolors=[colors[0],'none'],alpha=0.5) # Least extreme 0.003 & 0.01
        axes11.scatter([0,1],[S8s_pca_errors[1],S8s_pca_errors[4]]        ,color=[colors[1],colors[1]],marker=fmts[1],facecolors=[colors[1],'none'],alpha=0.5) # Intermediate  0.003 & 0.01
        axes11.scatter([0,1],[S8s_pca_errors[2],S8s_pca_errors[5]]        ,color=[colors[2],colors[2]],marker=fmts[2],facecolors=[colors[2],'none'],alpha=0.5) # Most extreme  0.003 & 0.01

        # Axes labels 
        axes00.set_ylabel(r"$\Omega_m$ bias, error")
        axes10.set_ylabel(r"$S_8$ bias, error")
        axes10.set_xlabel(r"# of PCs")
        axes10.set_xticks([0, 1, 2, 3, 4, 5])
        axes10.set_xlim(-0.2,5.2)
        axes01.set_title(r"Shift $\Delta_z^i$")
        axes11.set_xlabel(r"Prior")
        axes11.set_xticks([0, 1])
        axes11.set_xticklabels(['a', 'b'])
        axes00.xaxis.set_visible(False)
        for ax in [axes01]:
            ax.yaxis.set_visible(False)
            ax.xaxis.set_visible(False)
        for ax in [axes11]:
            ax.yaxis.set_visible(False)

        # Legend for marker types
        from matplotlib.lines import Line2D
        marker_labels = {'o': 'Least Extreme', '^': 'Intermediate', 's': 'Most Extreme'}
        colors = ["#1b5f6f", "#E69F00", "#882255"]

        legend_handles = [
        Line2D([0], [0], marker=m, color=c, linestyle='None', markersize=8, label=l)
        for (m, l), c in zip(marker_labels.items(), colors)]

        # Add line-style legend entries
        line_styles = {'-': 'Bias','--': r'$1\sigma$ Error',}
        line_handles = [
        Line2D([0], [0], color='gray', linestyle=ls, linewidth=2, label=lbl)
        for ls, lbl in line_styles.items()]

        # Combine both lists
        legend_handles.extend(line_handles)

        # Place legend above axes[0]
        axes00.legend(
        handles=legend_handles,
        ncol=5,               # horizontal layout
        title=None,           # no title
        fontsize=10.8,
        loc='lower center',   # anchor at bottom center of bbox_to_anchor
        bbox_to_anchor=(0.5, 1.0)  # x=0.5 centered, y=1.05 slightly above axes
        )

        fig_name="bias_error_mitigation_least_intermediate_most_extreme.pdf"
        plt.savefig(fig_name)  
        print(fig_name)

    return None

def get_inline_latex_w0wa():
    print(f"EXECUTING get_inline_latex_w0wa @ LINE {inspect.currentframe().f_lineno}")
    NPCs = 5
    combined_indx = list(range(0,NPCs+1)) + \
                    list(range(21,21+NPCs+1)) + \
                    list(range(42,42+NPCs+1))

    chains = [] 
    chains_baseline = [] 
    label = ''
    pattern = r"=\s*(-?[\d.]+)\^\{\+([\d.]+)\}_\{\-([\d.]+)\}"

    fig = plt.figure(figsize=(10, 4), constrained_layout=True)
    gs = fig.add_gridspec(4, 1, wspace=0.0)

    axes0 = fig.add_subplot(gs[0])
    axes1 = fig.add_subplot(gs[1],sharex=axes0) 
    axes2 = fig.add_subplot(gs[2],sharex=axes1)
    axes3 = fig.add_subplot(gs[3],sharex=axes2) 

    for i in combined_indx:
        samples = loadMCSamples(f"{chains_path}/weighted_pca_w0wa_Model_sc1bd4_Fiducial_sc1bd4/MCMC{i}", settings=settings)
        sigma8,omegam = samples["sigma8"],samples["omegam"]
        S8 = sigma8*np.sqrt(omegam/0.3)
        samples.addDerived(S8,name="S8",label="S_8")
        chains.append(samples) 

    omegams_pca = []
    S8s_pca = []
    w0s = []
    was = []
    # PCs
    for ch in chains:
        omegam_stats = ch.getInlineLatex("omegam", limit=2)
        S8_stats = ch.getInlineLatex("S8", limit=2)
        w0s_stats = ch.getInlineLatex("w", limit=2)
        was_stats = ch.getInlineLatex("wa", limit=2)
        print("PCs:",omegam_stats,S8_stats,w0s_stats,was_stats)                    
        match0 = re.search(pattern, omegam_stats)
        match1 = re.search(pattern, S8_stats)
        match2 = re.search(pattern, w0s_stats)
        match3 = re.search(pattern, was_stats)
        if match0:
            central = float(match0.group(1))
            upper = float(match0.group(2))
            lower = -float(match0.group(3))  # make it negative explicitly
            omegams_pca.append([central,lower,upper])
        if match1:
            central = float(match1.group(1))
            upper = float(match1.group(2))
            lower = -float(match1.group(3))  # make it negative explicitly
            S8s_pca.append([central,lower,upper])
        if match2:
            central = float(match2.group(1))
            upper = float(match2.group(2))
            lower = -float(match2.group(3))  # make it negative explicitly
            w0s.append([central,lower,upper])
        if match3:
            central = float(match3.group(1))
            upper = float(match3.group(2))
            lower = -float(match3.group(3))  # make it negative explicitly
            was.append([central,lower,upper])
    
    omegam_fid = 0.3
    S8_fid = 0.8277*(0.3/0.3)**0.5
    w0_fid = -1
    wa_fid = 0

    axes0.axhline(y=omegam_fid, color='gray', linestyle='--', linewidth=1.2)
    axes1.axhline(y=S8_fid    , color='gray', linestyle='--', linewidth=1.2)
    axes2.axhline(y=w0_fid    , color='gray', linestyle='--', linewidth=1.2)
    axes3.axhline(y=wa_fid    , color='gray', linestyle='--', linewidth=1.2)

    ecolors = ["#1b5f6f"]*(NPCs+1)+\
              ["#E69F00"]*(NPCs+1)+\
              ["#882255"]*(NPCs+1)
    
    fmts = ["o"]*(NPCs+1)+\
           ["^"]*(NPCs+1)+\
           ["s"]*(NPCs+1)
    
    dx = [-.1]*(NPCs+1)+\
            [0]*(NPCs+1)+\
            [.1]*(NPCs+1)

    cases = list(range(len(combined_indx)))
    npcs = list(range(NPCs+1))
    pairs = [(npcs[i % len(npcs)], c) for i, c in enumerate(cases)]
    print("pairs:",pairs)
    print("omegam:",omegams_pca)
    print("S8:",S8s_pca)
    print("w0:",w0s)
    print("wa:",was)
    for p in pairs:
        npc = p[0]
        case = p[1]
        axes0.errorbar(x=npc+dx[case], y=omegams_pca[case][0], yerr=[[abs(omegams_pca[case][1])], [omegams_pca[case][2]]],
                        fmt=fmts[case], capsize=4, color=ecolors[case],ecolor=ecolors[case], elinewidth=1.2)
        axes1.errorbar(x=npc+dx[case], y=S8s_pca[case][0], yerr=[[abs(S8s_pca[case][1])], [S8s_pca[case][2]]],
                        fmt=fmts[case], capsize=4, color=ecolors[case], ecolor=ecolors[case], elinewidth=1.2)
        axes2.errorbar(x=npc+dx[case], y=w0s[case][0], yerr=[[abs(w0s[case][1])], [w0s[case][2]]],
                        fmt=fmts[case], capsize=4, color=ecolors[case], ecolor=ecolors[case], elinewidth=1.2)
        axes3.errorbar(x=npc+dx[case], y=was[case][0], yerr=[[abs(was[case][1])], [was[case][2]]],
                        fmt=fmts[case], capsize=4, color=ecolors[case], ecolor=ecolors[case], elinewidth=1.2)
    
    ecolors = ["#1b5f6f","#E69F00","#882255"]
    fmts = ["o","^","s"]
    dx = [-.01,0,.01]
    
    print(omegams_pca)
    print(S8s_pca)
    print(w0s)
    print(was)
    axes0.set_ylabel(r"$\Omega_m$")
    axes1.set_ylabel(r"$S_8$")
    axes2.set_ylabel(r"$w_0$")
    axes3.set_ylabel(r"$w_a$")
    axes3.set_xlabel(r"# of PCs")
    axes3.set_xticks([0, 1, 2, 3, 4, 5])
    axes3.set_xlim(-0.2,5.2)
    axes0.xaxis.set_visible(False)
    axes1.xaxis.set_visible(False)
    axes2.xaxis.set_visible(False)
    # for ax in [axes01, axes11]:
    #     ax.yaxis.set_visible(False)
    #     ax.xaxis.set_visible(False)    

    # Legend for marker types
    from matplotlib.lines import Line2D
    marker_labels = {'o': 'Least Extreme', '^': 'Intermediate', 's': 'Most Extreme'}
    colors = ["#1b5f6f", "#E69F00", "#882255"]

    legend_handles = [
    Line2D([0], [0], marker=m, color=c, linestyle='None', markersize=8, label=l)
    for (m, l), c in zip(marker_labels.items(), colors)]

    # Place legend above axes[0]
    axes0.legend(
    handles=legend_handles,
    ncol=3,               # horizontal layout
    title=None,           # no title
    fontsize=13,
    loc='lower center',   # anchor at bottom center of bbox_to_anchor
    bbox_to_anchor=(0.5, 1.0)  # x=0.5 centered, y=1.05 slightly above axes
    )
    fig_name="mitigation_least_intermediate_most_extreme_w0wa.pdf"
    plt.savefig(fig_name) # Paper Figure 12 
    print(fig_name)
    return None

def triangle_alphas_non_mixing():
    print(f"EXECUTING triangle_alphas_non_mixing @ LINE {inspect.currentframe().f_lineno}")
    g = plots.get_subplot_plotter()
    g.settings.title_limit_fontsize = 13
    g.settings.legend_fontsize = 20
    g.settings.axes_labelsize = 20
    g.settings.axes_fontsize = 15    
    chains = [] #initial_chains.copy()
    label=''
    bconfig = {"legend_ncol":1,
               "params":["roman_alpha_1",
                         "roman_alpha_2",
                         "roman_alpha_3",
                         "roman_alpha_4",
                         "roman_alpha_5"],
                "legend_labels":["Least Extreme",
                                 "Intermediate",
                                 "Most Extreme"], 
                "markers":{"roman_alpha_1": 0,
                           "roman_alpha_2": 0,
                           "roman_alpha_3": 0,
                           "roman_alpha_4": 0,
                           "roman_alpha_5": 0,
                           },                         
                "line_args":[{"ls": "-" , "color": "#1b5f6f", "lw":1.5},
                             {"ls": "--", "color": "#E69F00", "lw":1.5},
                             {"ls": "-.", "color": "#882255", "lw":1.5},],
                "contour_colors":["#1b5f6f","#E69F00","#882255"],
                "contour_ls":["-","--","-."],
                "contour_lws":[1.5,1.5,1.5],
                "filled":False}        
    for i in [86,107,128]:
        samples = loadMCSamples(f"{chains_w_path}{i}", settings=settings)
        chains.append(samples) 
        label+=f"{i},"
    g.triangle_plot(roots=chains,**bconfig)
    fig_name=f"triangle_plot_weighted_alphas_{label}_Model_sc1bd4_Fiducial_sc1bd4.pdf"
    g.export(fig_name)
    print(fig_name)
    return None

################################################
### MIXING SCENARIOS - PAPER SECTION "RESULTS"
################################################

def triangle_plot_mixing(plot_type="Omegam_x_S8",mod="sc1bd5",fid="sc1bd4"):
    print(f"EXECUTING triangle_plot_mixing @ LINE {inspect.currentframe().f_lineno}")
    g = plots.get_subplot_plotter()
    g.settings.title_limit_fontsize = 13
    g.settings.legend_fontsize = 20
    g.settings.axes_labelsize = 20
    g.settings.axes_fontsize = 15
    colors1 = [
        "#1b5f6f",  # dark teal
        "#E69F00",  # orange
        "#882255",  # wine
        "#56B4E9",  # sky blue
        "#009E73",  # bluish green
        "black"  ,  # black
        "#0072B2",  # blue
        "#D55E00",  # vermillion
        "#CC79A7",  # reddish purple
        "#999999",  # grey
        "#117733",  # dark green
    ]

    chains = [] #initial_chains.copy()
    label=''
    if plot_type=="Omegam_x_S8":
        g.settings.legend_fontsize = 12
        bconfig = {"params":["omegam","S8"],
                   "markers":{"omegam": 0.3, "S8":0.8277*(0.3/0.3)**0.5},
                   "legend_ncol":1,
                   "legend_labels":["0 PC","1 PC","2 PC","3 PC","4 PC","5 PC"],
                   "line_args":[{"ls": "-" , "color": colors1[0], "lw":1.5},
                                {"ls": "--", "color": colors1[1], "lw":1.5},
                                {"ls": "-.", "color": colors1[2], "lw":1.5},
                                {"ls": ":" , "color": colors1[3], "lw":1.5},
                                {"ls": "-" , "color": colors1[4], "lw":1.5},
                                {"ls": "--", "color": colors1[5], "lw":1.5},
                                ],
                   "contour_colors":[colors1[0],colors1[1],colors1[2],colors1[3],colors1[4],colors1[5]],
                   "contour_ls":["-","--","-.",":","-","--"],
                   "contour_lws":[1.5,1.5,1.5,1.5,1.5,1.5],
                   "filled":False}        
        NPCs = 5
        for i in list(range(NPCs+1)):
            samples = loadMCSamples(f"{chains_path}/weighted_pca_Model_{mod}_Fiducial_{fid}/MCMC{i}", settings=settings)
            sigma8,omegam = samples["sigma8"],samples["omegam"]
            S8 = sigma8*np.sqrt(omegam/0.3)
            samples.addDerived(S8,name="S8",label="S_8")
            chains.append(samples) 
            label+=f"{i},"
        g.triangle_plot(roots=chains,**bconfig)
        fig_name=f"triangle_plot_weighted_{label}_Model_{mod}_Fiducial_{fid}.pdf"
        g.export(fig_name)
        g.export(fig_name)
        print(fig_name)
        ####
    elif plot_type=="roman_alphas":
        bconfig = {"params":["roman_alpha_1",
                             "roman_alpha_2",
                             "roman_alpha_3",
                             "roman_alpha_4",
                             "roman_alpha_5"
                            ],
                    "markers":{"roman_alpha_1": 0,
                               "roman_alpha_2": 0,
                               "roman_alpha_3": 0,
                               "roman_alpha_4": 0,
                               "roman_alpha_5": 0,
                               },        
                   "legend_ncol":1,
                   "line_args":[
                            {"ls": "-"  , "color": colors1[0], "lw":1.5},
                            {"ls": "--" , "color": colors1[1], "lw":1.5},
                            {"ls": ":"  , "color": colors1[2], "lw":1.5},
                            ],
                   "contour_colors":[colors1[0],colors1[1],colors1[2]],
                   "contour_ls":["-","--",":"],
                   "contour_lws":[1.5,1.5,1.5],
                   "filled":False}        
        for i in [5]:
            for p in [chains_sc1bd5_x_sc1bd4_path,
                      chains_sc1bd6_x_sc1bd4_path,
                      chains_sc1bd7_x_sc1bd4_path]:
                samples = loadMCSamples(f"{p}{i}", settings=settings)
                chains.append(samples) 
            label+=f"{i},"
        g.triangle_plot(roots=chains,
                        **bconfig,
                        legend_labels=["PCs: DRM-D2\nFiducial datavector: DRM-D1",
                                       "PCs: DRM-D3\nFiducial datavector: DRM-D1",
                                       "PCs: DRM-D4\nFiducial datavector: DRM-D1",
                                       ]
                        )
        fig_name=f"triangle_plot_weighted_alphas_{label}_Model_{mod}_Fiducial_{fid}.pdf"
        g.export(fig_name)
        print(fig_name)
    return None

def get_inline_latex_mixing_scenarios():
    print(f"EXECUTING get_inline_latex_mixing_scenarios @ LINE {inspect.currentframe().f_lineno}")
    chains_model_x_fiducial_path = lambda model,fiducial: f'{chains_path}/weighted_pca_Model_{model}_Fiducial_{fiducial}/MCMC'
    scenarios = [("sc1bd5","sc1bd4"),("sc1bd6","sc1bd4"),("sc1bd7","sc1bd4"),
                 ("sc2bd5","sc2bd4"),("sc2bd6","sc2bd4"),("sc2bd7","sc2bd4"),
                 ("sc3bd5","sc3bd4"),
                 ("sc3bd7","sc3bd4")
                 ]

    NPCs = 5
    print("len(scenarios): ",len(scenarios))
    print("NPCs: ",NPCs)
    idxs = list(range(0,NPCs+1))

    pattern = r"=\s*([\d.]+)\^\{\+([\d.]+)\}_\{\-([\d.]+)\}"

    fig = plt.figure(figsize=(15, 4), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, wspace=0.0, hspace=0.0)

    axes0 = fig.add_subplot(gs[0])
    axes1 = fig.add_subplot(gs[1],sharex=axes0)

    colors1 = [
        "#1b5f6f",  # dark teal
        "#E69F00",  # orange
        "#882255",  # wine
        "#56B4E9",  # sky blue
        "#009E73",  # bluish green
        "#F0E442",  # yellow
        "#0072B2",  # blue
        "#D55E00",  # vermillion
        "#CC79A7",  # reddish purple
        "#117733",  # dark green
        "#999999",  # grey
    ]
    
    fmts = {"sc1bd4":"o" , "sc2bd4":"s" , "sc3bd4":"^"}

    ecolors = {
               "sc1bd5":colors1[0],"sc1bd6":colors1[1],"sc1bd7":colors1[2],
               "sc2bd5":colors1[3],"sc2bd6":colors1[4],"sc2bd7":colors1[5],
               "sc3bd5":colors1[6],"sc3bd7":colors1[7]
               }
    
    shifts = {
               "sc1bd5":-1.6,"sc1bd6":-1.5,"sc1bd7":-1.4,
               "sc2bd5":-0.8,"sc2bd6":-0.7,"sc2bd7":-0.6,
               "sc3bd5":-0.1,"sc3bd7":0
               }

    x_positions = {
                   # Fiducial: Model: N of PCs: ...
                   "sc1bd4":{"sc1bd5":{"npc":{"0":0.0,"1":2.6,"2":5.1,"3":7.5,"4":10.0,"5":12.5}},
                             "sc1bd6":{"npc":{"0":0.1,"1":2.7,"2":5.2,"3":7.6,"4":10.1,"5":12.6}},
                             "sc1bd7":{"npc":{"0":0.2,"1":2.8,"2":5.3,"3":7.7,"4":10.2,"5":12.7}},
                             },
                   "sc2bd4":{"sc2bd5":{"npc":{"0":0.7,"1":3.3,"2":5.7,"3":8.2,"4":10.7,"5":13.2}},
                             "sc2bd6":{"npc":{"0":0.8,"1":3.4,"2":5.8,"3":8.3,"4":10.8,"5":13.3}},
                             "sc2bd7":{"npc":{"0":0.9,"1":3.5,"2":5.9,"3":8.4,"4":10.9,"5":13.4}},
                             },
                   "sc3bd4":{"sc3bd5":{"npc":{"0":1.5,"1":4.0,"2":6.4,"3":8.9,"4":11.4,"5":13.9}},
                             "sc3bd7":{"npc":{"0":1.6,"1":4.1,"2":6.5,"3":9.0,"4":11.5,"5":14.0}},
                             },
                   }

    chains = [] 
    info_chains = []
    for mod,fid in scenarios:
        print(mod,fid)
        for i in idxs:
            samples = loadMCSamples(f"{chains_model_x_fiducial_path(mod,fid)}{i}", settings=settings)
            sigma8,omegam = samples["sigma8"],samples["omegam"]
            S8 = sigma8*np.sqrt(omegam/0.3)
            samples.addDerived(S8,name="S8",label="S_8")
            chains.append(samples) 
            info_chains.append((samples,i,mod,fid,ecolors[mod],fmts[fid],shifts[mod]))

    omegams_pca = []
    S8s_pca = []
    w0s = []
    was = []
    for chi,ch in enumerate(chains):
        omegam_stats = ch.getInlineLatex("omegam", limit=2)
        S8_stats = ch.getInlineLatex("S8", limit=2)
        match1 = re.search(pattern, omegam_stats)
        match2 = re.search(pattern, S8_stats)
        print("Chain index, Statistics:",chi,omegam_stats,S8_stats)
        if match1:
            central = float(match1.group(1))
            upper = float(match1.group(2))
            lower = -float(match1.group(3))  # make it negative explicitly
            omegams_pca.append([central,lower,upper])
        if match2:
            central = float(match2.group(1))
            upper = float(match2.group(2))
            lower = -float(match2.group(3))  # make it negative explicitly
            S8s_pca.append([central,lower,upper])               
    
    print("omegams_pca: ",omegams_pca)
    print("S8s_pca: ",S8s_pca)

    omegam_fid = 0.3
    S8_fid = 0.8277*(0.3/0.3)**0.5

    axes0.axhline(y=omegam_fid, color='gray', linestyle='--', linewidth=1.2)
    axes1.axhline(y=S8_fid    , color='gray', linestyle='--', linewidth=1.2)
    shift_global = 0
    print("info_chains: ",info_chains)
    for i,info_chain in enumerate(info_chains):
        samples,npc,mod,fid,ecolor,fmt,shift = info_chain
        x_pos = x_positions[fid][mod]["npc"][f"{npc}"]
        axes0.errorbar(x=x_pos, y=omegams_pca[i][0], yerr=[[abs(omegams_pca[i][1])], [omegams_pca[i][2]]],fmt=fmt, color=ecolor, ecolor=ecolor)
        axes1.errorbar(x=x_pos, y=S8s_pca[i][0]    , yerr=[[abs(S8s_pca[i][1])]    , [S8s_pca[i][2]]]    ,fmt=fmt, color=ecolor, ecolor=ecolor)
    
    axes0.set_ylabel(r"$\Omega_m$")
    axes1.set_ylabel(r"$S_8$")
    axes1.set_xlabel(r"# of PCs")
    x_lim_a = x_positions["sc1bd4"]["sc1bd5"]["npc"]["0"] # first point
    x_lim_b = x_positions["sc3bd4"]["sc3bd7"]["npc"]["5"] # last point
    axes1.set_xlim(x_lim_a-0.1,x_lim_b+0.1)
    axes0.xaxis.set_visible(False)

    tick_positions = np.linspace(-0.1, 14, 6)
    fid_ref, mod_ref = "sc2bd4", "sc2bd6"
    tick_positions = [x_positions[fid_ref][mod_ref]["npc"][f"{i}"] for i in range(NPCs+1)]
    tick_labels = [0, 1, 2, 3, 4, 5]

    axes1.set_xticks(tick_positions)
    axes1.set_xticklabels(tick_labels)

    # Legend for marker types
    from matplotlib.lines import Line2D
    marker_labels = {'o': 'M1-D1', 's': 'M2-D1', '^': 'M3-D1'}
    colors = ["gray"]*3

    legend_handles1 = [
    Line2D([0], [0], marker=m, color=c, linestyle='None', markersize=8, label=l)
    for (m, l), c in zip(marker_labels.items(), colors)]

    # Place legend above axes[0]
    leg1 = axes0.legend(
        handles=legend_handles1,
        ncol=3,
        fontsize=12,
        loc='lower left',
        bbox_to_anchor=(0, 1.0),
    )

    # --- Legend 2: color meanings (for example) ---
    color_labels = {colors1[0]: 'M1-D2',colors1[1]: 'M1-D3',colors1[2]: 'M1-D4',
                    colors1[3]: 'M2-D2',colors1[4]: 'M2-D3',colors1[5]: 'M2-D4',
                    colors1[6]: 'M3-D2',colors1[7]: 'M3-D4'}
    legend_handles2 = [
        Line2D([0], [0], color=c, linestyle='-', label=l, lw=2)
        for c, l in color_labels.items()]

    # Add second legend to the same axes
    leg2 = axes0.legend(
        handles=legend_handles2,
        ncol=8,
        fontsize=12,
        loc='lower right',
        bbox_to_anchor=(1.0, 1.0),
    )

    axes0.add_artist(leg1)

    fig_name="mitigation_mixing_scenarios.pdf"
    plt.savefig(fig_name) 
    print(fig_name)
    return None


################################
###### FUNCTION INVENTORY ######
################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--func", required=True, help="Function name to run")
    args = parser.parse_args()

    # Map function names to actual functions
    func_map = {
        "triangle_plot":triangle_plot,
        "plot_chisq":plot_chisq,
        "plot_2d":plot_2d,
        "get_inline_latex_lcdm":get_inline_latex_lcdm,
        "bias_error_mitigation_lcdm":bias_error_mitigation_lcdm,
        "get_inline_latex_w0wa":get_inline_latex_w0wa,
        "triangle_alphas_non_mixing":triangle_alphas_non_mixing,
        # "triangle_plot_mixing":triangle_plot_mixinglot_type="roman_alphas",combination="Model_sc1bd6_Fiducial_sc1bd4"),
        # "triangle_plot_mixing":triangle_plot_mixinglot_type="Omegam_x_S8",mod="sc3bd7",fid="sc3bd4"),
        "triangle_plot":triangle_plot,
        "get_inline_latex_mixing_scenarios":get_inline_latex_mixing_scenarios,
    }

    if args.func not in func_map:
        raise ValueError(f"Unknown function '{args.func}'. Available: {list(func_map)}")

    # Call the selected function
    func_map[args.func]()