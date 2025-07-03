from getdist import plots,loadMCSamples,MCSamples

burnin=0.5
settings = {"ignore_rows": burnin}
path = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/results/chains/roman_sc1bd4_g/"

no_pc = loadMCSamples(path+"NZ_MCMC0", settings=settings)
one_pc = loadMCSamples(path+"NZ_MCMC1", settings=settings)
two_pc = loadMCSamples(path+"NZ_MCMC2", settings=settings)
three_pc = loadMCSamples(path+"NZ_MCMC3", settings=settings)

# chain.getInlineLatex("x1", limit=1)

g = plots.get_subplot_plotter()
g.plots_1d([no_pc,one_pc,two_pc,three_pc],["omegam","sigma8","w","wa"],
           markers={"omegam":0.3,"sigma8":0.8120,"w":-1,"wa":0}, # planck best-fit 1807.06209
           nx=4,legend_ncol=4,legend_labels=['No PC', '1 PC', '2 PC', '3 PC'],
           colors=['#1b5f6f',"#E69F00","#56B4E9","#D55E00"],
           ls=['-','--',':','-.'],lws=[3,3,3,3])
g.export("test.pdf")