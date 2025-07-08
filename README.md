# cocoa_photoz  
MCMC PCA SCHEDULE PROGRAM  

# PCA Method  
A: PCs from Cn  
B: PCs from nz_compression  

# Dv used  
same: dv generated with same params as used in CosmoCov  
diff: dv generated with diff params as used in CosmoCov  

# Prior on alpha:  
Gauss: Gaussian (0,1)  
Flat: Flat (-3,3)  

#### MCMC  
### Prior: Gauss  
## Method A: PCs from Cn  
# dv: same  
0 / running  
1 / running  
2 / running  
3 / running  
# dv: diff  
4 / running  
5 / running  
6 / running  
7 / running  
## Method B: PCs from nz_compression  
# dv: same  
8 = 0  
9 / running  
10 / running  
11 / running  
# dv: diff  
12 = 4  
13 / running  
14 / running  
15 / running  
### Prior: Flat  