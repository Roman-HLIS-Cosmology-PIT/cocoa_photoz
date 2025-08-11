Steps to run an evaluation and a mcmc with the PCA model for n(z)  
Assuming you already have CoCoA installed, follow these steps:  

1) cd cocoa/Cocoa  
2) git clone git@github.com:diogohf/cocoa_photoz.git  
3) Copy and past the following files into the Roman Real folders likelihood, interface, and data:  
    cp cocoa_photoz/modifications/roman_real/likelihood/* projects/roman_real/likelihood/  
    cp cocoa_photoz/modifications/roman_real/interface/* projects/roman_real/interface/  
    cp -r cocoa_photoz/modifications/roman_real/data/* projects/roman_real/data/  
4) Decompress the covarianc matrix: gunzip projects/roman_real/data/sc1bd4/cov_sc1bd4r0.gz
5) Compile Roman Real: source projects/roman_real/scripts/compile_roman_real.sh   


4) Minimal MCMC example with PCA:  
MCMC index / Roman scenario / realization index / chain status (NS=Not Started, R=Running, C=Converged)    
### Setup: No PC & No Shift. Synthetic data vector: ξ±(nz_fid = nz_991213). Model: ξ±(nz = nz_fid)  
57         / sc1bd4         / 991213            / NS        
### Setup: No PC & Shift. Synthetic data vector: ξ±(nz_fid = nz_991213). Model: ξ±(nz = nz_bar). nz_bar is the mean n(z) of 10^6 realization of the scenario sc1bd4  
59         / sc1bd4         / 991213            / NS        
### Setup: 1 PC & No Shift. Synthetic data vector: ξ±(nz_fid = nz_991213). Model: ξ±(nz = nz_bar + α_1*PC_1).
61         / sc1bd4         / 991213            / NS        
### Setup: 5 PC & No Shift. Synthetic data vector: ξ±(nz_fid = nz_991213). Model: ξ±(nz = nz_bar + α_1 * PC_1 +...+ α_5 * PC_5).
69         / sc1bd4         / 991213            / NS        
### Setup: 10 PC & No Shift. Synthetic data vector: ξ±(nz_fid = nz_991213). Model: ξ±(nz = nz_bar + α_1 * PC_1 +...+ α_10 * PC_10).
79         / sc1bd4         / 991213            / NS        