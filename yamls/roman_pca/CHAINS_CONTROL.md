# cocoa_photoz  
MCMC PCA SCHEDULE PROGRAM   

# Model  
lcdm  
wcdm   
w0wa      

# Prior on alpha:  
Gaussian (0,1)  

# ###################################################################  
# ###################################################################  

# UNWEIGHTED PCA LCDM 

## TEST 1 - fiducial cosmology from planck 18  
0  / sc1bd4 + nbar / no pc + shift     / running  
1  / sc1bd4 + nbar / no pc & no shift  / C  
2  / sc1bd4 + nbar / 1 pc              / C    
3  / sc1bd4 + nbar / 2 pc              / C   
4  / sc1bd4 + nbar / 3 pc              / C   
5  / sc1bd4 + nbar / 4 pc              / C   
6  / sc1bd4 + nbar / 5 pc              / C   
7  / sc1bd4 + nbar / 6 pc              / running   
8  / sc1bd4 + nbar / 7 pc              / C   
9  / sc1bd4 + nbar / 8 pc              / C   
10 / sc1bd4 + nbar / 9 pc              / running   
11 / sc1bd4 + nbar / 10 pc             / C  
12 / sc1bd4 + nzr0 / no pc + shift     /  *rerun*  
13 / sc1bd4 + nzr0 / no pc & no shift  /  running  

## TEST 2 - fiducial cosmology from planck 18  
14 sc1bd4 + nzr0 / no pc + shift / same as 12  
15 sc1bd5 + nzr0 / no pc + shift / running on edge  
16 sc1bd6 + nzr0 / no pc + shift / running on edge  
17 sc1bd7 + nzr0 / no pc + shift / running on edge  
18 sc2bd4 + nzr0 / no pc + shift / running on edge  
19 sc2bd5 + nzr0 / no pc + shift / running on edge  
20 sc2bd6 + nzr0 / no pc + shift / running on edge  
21 sc2bd7 + nzr0 / no pc + shift / running on edge  
22 sc3bd4 + nzr0 / no pc + shift / running on edge  
23 sc3bd5 + nzr0 / no pc + shift / running on edge  
24 sc3bd7 + nzr0 / no pc + shift / running on edge  
## TEST 3 - fiducial cosmology from roman real  
25 sc1bd4 / !file not created yet!   
26 sc1bd5 / !file not created yet!   
27 sc1bd6 / !file not created yet!   
28 sc1bd7 / !file not created yet!   
29 sc2bd4 / !file not created yet!   
30 sc2bd5 / !file not created yet!   
31 sc2bd6 / !file not created yet!   
32 sc2bd7 / !file not created yet!   
33 sc3bd4 / !file not created yet!   
34 sc3bd5 / !file not created yet!   
35 sc3bd7 / !file not created yet!   

## TEST 4 - fiducial cosmology from roman real: see l2_norm_for_nz.py  
### No PC + Shift with synthetic xi(nfid) + model xi(nfid)  
36 sc1bd4 [extreme realization = 997946] / PAUSED (jpl)        [scp to edge done]  
37 sc1bd4 [extreme realization = 870731] / canceled (jpl)       [scp to edge done]    
38 sc1bd4 [extreme realization = 174820] / canceled (jpl)       [scp to edge done]    
39 sc1bd4 [extreme realization = 991213] / PAUSED (jpl)        [scp to edge done]    
40 sc1bd4 [extreme realization = 645431] / canceled (jpl)       [scp to edge done]    
41 sc1bd4 [extreme realization = 131103] / canceled (jpl)       [scp to edge done]    
42 sc1bd4 [extreme realization = 808898] / canceled (edge)      [scp to jpl done]      
43 sc1bd4 [extreme realization = 682038] / canceled (edge)      [scp to jpl done]      
44 sc1bd4 [extreme realization = 472580] / canceled (jpl)       [scp to edge done]    
45 sc1bd4 [extreme realization = 92550 ] / canceled (jpl)       [scp to edge done]     

### No PC + Shift with synthetic xi(nfid) + model xi(nbar)
46 sc1bd4 [realization = 997946] / PAUSED (edge)       
47 sc1bd4 [realization = 870731] / canceled (edge)      [scp to jpl done]   
48 sc1bd4 [realization = 174820] / canceled (edge)      [scp to jpl done]     
49 sc1bd4 [realization = 991213] / PAUSED (edge)         
50 sc1bd4 [realization = 645431] / canceled (edge)      [scp to jpl done]     
51 sc1bd4 [realization = 131103] / canceled (edge)      [scp to jpl done]     
52 sc1bd4 [realization = 808898] / canceled (edge)      [scp to jpl done]     
53 sc1bd4 [realization = 682038] / canceled (edge)      [scp to jpl done]     
54 sc1bd4 [realization = 472580] / canceled (edge)      [scp to jpl done]     
55 sc1bd4 [realization = 92550 ] / canceled (edge)      [scp to jpl done]     

### No PC + No Shift with synthetic xi(nfid) + model xi(fid)
56 sc1bd4 [realization = 997946] / C (jpl)     
57 sc1bd4 [realization = 991213] / C (jpl)     
### x PC + No Shift with synthetic xi(nfid) + model xi(nbar)
58 / 0 PC / sc1bd4 [realization = 997946] / C (jpl)     
59 / 0 PC / sc1bd4 [realization = 991213] / C (jpl)     
60 / 1 PC / sc1bd4 [realization = 997946] / C (jpl)     
61 / 1 PC / sc1bd4 [realization = 991213] / C (edge)     
62 / 2 PC / sc1bd4 [realization = 997946] / RUNNING (edge)     
63 / 2 PC / sc1bd4 [realization = 991213] / C (edge)     
64 / 3 PC / sc1bd4 [realization = 997946] / C (edge)     
65 / 3 PC / sc1bd4 [realization = 991213] / C (edge)     
66 / 4 PC / sc1bd4 [realization = 997946] / RUNNING (edge)     
67 / 4 PC / sc1bd4 [realization = 991213] / C (edge)     
68 / 5 PC / sc1bd4 [realization = 997946] / RUNNING (edge)     
69 / 5 PC / sc1bd4 [realization = 991213] / C (edge)     
70 / 6 PC / sc1bd4 [realization = 997946] / RUNNING (edge)     
71 / 6 PC / sc1bd4 [realization = 991213] / C (jpl)     
72  / 7 PC / sc1bd4 [realization = 997946] / RUNNING (jpl)     
73  / 7 PC / sc1bd4 [realization = 991213] / C (jpl)     
74 / 8 PC / sc1bd4 [realization = 997946] / WAITING (jpl)     
75 / 8 PC / sc1bd4 [realization = 991213] / C (jpl)     
76 / 9 PC / sc1bd4 [realization = 997946] / waiting (jpl)     
77 / 9 PC / sc1bd4 [realization = 991213] / RUNNING (jpl)     
78 / 10 PC / sc1bd4 [realization = 997946] / RUNNING (jpl)     
79 / 10 PC / sc1bd4 [realization = 991213] / C (jpl)       
80 / 11 PC / sc1bd4 [realization = 997946] / PAUSED (edge)     
81 / 11 PC / sc1bd4 [realization = 991213] / PAUSED (edge)      
82 / 12 PC / sc1bd4 [realization = 997946] / PAUSED (edge)     
83 / 12 PC / sc1bd4 [realization = 991213] / WAITING (edge)      
84 / 13 PC / sc1bd4 [realization = 997946] / PAUSED (edge)     
85 / 13 PC / sc1bd4 [realization = 991213] / PAUSED (edge)      
86 / 14 PC / sc1bd4 [realization = 997946] / PAUSED (edge)     
87 / 14 PC / sc1bd4 [realization = 991213] / PAUSED (edge)      
88 / 15 PC / sc1bd4 [realization = 997946] / PAUSED (edge)           
89 / 15 PC / sc1bd4 [realization = 991213] / PAUSED (edge)      
90 / 16 PC / sc1bd4 [realization = 997946] / PAUSED (edge)      
91 / 16 PC / sc1bd4 [realization = 991213] / PAUSED (edge)      
92 / 17 PC / sc1bd4 [realization = 997946] / PAUSED (edge)     
93 / 17 PC / sc1bd4 [realization = 991213] / PAUSED (edge)      
94 / 18 PC / sc1bd4 [realization = 997946] / PAUSED (edge)     
95 / 18 PC / sc1bd4 [realization = 991213] / PAUSED (edge)      
96 / 19 PC / sc1bd4 [realization = 997946] / PAUSED (edge)     
97 / 19 PC / sc1bd4 [realization = 991213] / PAUSED (edge)      
98 / 20 PC / sc1bd4 [realization = 997946] / PAUSED (edge)     
99 / 20 PC / sc1bd4 [realization = 991213] / PAUSED (edge)      

# ###################################################################  
# ###################################################################  

# WEIGHTED PCA $\Lambda\mathrm{CDM}$

## PCs $\text{\textcolor{red}{trained}}$ on M1-D1 (sc1bd4), $\text{\textcolor{cyan}{reconstruct}}$ nz from M1-D1 (sc1bd4)  
(same as TEST 4 of unweighted PCA)  
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca/`  

$
\begin{align}
\boxed{\xi_\mathrm{model}=\xi_\mathrm{cocoa}[\bar{n}_\mathrm{DRM-D1}+\sum\alpha_i\mathrm{PC}_\mathrm{DRM-D1}]};
\end{align}
$     
With shift $\Delta_{z,i}$. Sinthetic data: $\boxed{\xi_\mathrm{fid}=\xi_\mathrm{cocoa}[n_{\mathrm{fid},\mathrm{DRM-D1}}^r]~;~r=991213}$  
MCMC_BASELINEi / number of Δs per bin / realization index / status    
external_nz_modeling = 0  
1 / 1 / [realization = 991213] / running (edge)   
2 / 2 / [realization = 991213] / running (edge)     
3 / 3 / [realization = 991213] / running (edge)     
4 / 4 / [realization = 991213] / running (edge)     
5 / 5 / [realization = 991213] / running (edge)     
6 / 6 / [realization = 991213] / running (edge)     
7 / 7 / [realization = 991213] / running (edge)     
8 / 8 / [realization = 991213] / running (edge)     
9 / 9 / [realization = 991213] / running (edge)      
external_nz_modeling = 1   
10 / 1 / [realization = 991213] / running (edge)  
11 / 2 / [realization = 991213] / running (edge)  
12 / 3 / [realization = 991213] / running (edge)  
13 / 4 / [realization = 991213] / running (edge)  
14 / 5 / [realization = 991213] / running (edge)  
15 / 6 / [realization = 991213] / running (edge)  
16 / 7 / [realization = 991213] / running (edge)  
17 / 8 / [realization = 991213] / running (edge)  
18 / 9 / [realization = 991213] / running (edge)   
external_nz_modeling = 1   
19 / 9 / [realization = 9085] (nbar_sc1bd4.nz) / running (edge) - same as 21, but with tighter prior on the ∆zs    
20 / 9 / [realization = 9085] (sc1bd4_r9085.nz) / running (jpl)      
21 / 9 / [realization = 9085] (nbar_sc1bd4.nz) / running (edge) - same as 19, but with larger prior on the ∆zs     

No shift ($\Delta_{z,i}=0$). Sinthetic data: $\boxed{\xi_\mathrm{fid}=\xi_\mathrm{cocoa}[n_{\mathrm{fid},\mathrm{DRM-D1}}^r]~;~r=997946~\text{or}~r=991213}$  

$\alpha\mathrm{s}\sim\mathrm{G}(0,1)$ [uninformative prior]  

0 / 1 PC / 60 sc1bd4 [realization = 997946] / RUNNING (edge)     
1 / 1 PC / 61 sc1bd4 [realization = 991213] / C (edge)     
2 / 2 PC / 62 sc1bd4 [realization = 997946] / C (edge)     
3 / 2 PC / 63 sc1bd4 [realization = 991213] / C (edge)     
4 / 3 PC / 64 sc1bd4 [realization = 997946] / C (edge)     
5 / 3 PC / 65 sc1bd4 [realization = 991213] / C (edge)     
6 / 4 PC / 66 sc1bd4 [realization = 997946] / C (edge)     
7 / 4 PC / 67 sc1bd4 [realization = 991213] / C (edge)     
8 / 5 PC / 68 sc1bd4 [realization = 997946] / C (edge)     
9 / 5 PC / 69 sc1bd4 [realization = 991213] / C (edge)           
10 / 6 PC / 70 sc1bd4 [realization = 997946] / C (edge)     
11 / 6 PC / 71 sc1bd4 [realization = 991213] / C (edge)     
12 / 7 PC / 72 sc1bd4 [realization = 997946] / C (edge)     
13 / 7 PC / 73 sc1bd4 [realization = 991213] / C (edge)     
14 / 8 PC / 74 sc1bd4 [realization = 997946] / C (edge)     
15 / 8 PC / 75 sc1bd4 [realization = 991213] / C (edge)     
16 / 9 PC / 76 sc1bd4 [realization = 997946] / C (edge)     
17 / 9 PC / 77 sc1bd4 [realization = 991213] / C (edge)     
18 / 10 PC / 78 sc1bd4 [realization = 997946] / C (jpl)            
19 / 10 PC / 79 sc1bd4 [realization = 991213] / C (jpl)                    
20 / 11 PC / 80 sc1bd4 [realization = 997946] / RUNNING (jpl)            
21 / 11 PC / 81 sc1bd4 [realization = 991213] / C (jpl)            
22/ 12 PC / 82 sc1bd4 [realization = 997946] / RUNNING (jpl)            
23/ 12 PC / 83 sc1bd4 [realization = 991213] / C (jpl)            
24 / 13 PC / 84 sc1bd4 [realization = 997946] / RUNNING (jpl)            
25 / 13 PC / 85 sc1bd4 [realization = 991213] / RUNNING (jpl)            
26 / 14 PC / 86 sc1bd4 [realization = 997946] / C (jpl)            
27 / 14 PC / 87 sc1bd4 [realization = 991213] / C (jpl)            
28 / 15 PC / 88 sc1bd4 [realization = 997946] / C (jpl)            
29 / 15 PC / 89 sc1bd4 [realization = 991213] / C (jpl)            
  
30 / 16 PC / 90 sc1bd4 [realization = 997946] / C (jpl)                   
31 / 16 PC / 91 sc1bd4 [realization = 991213] / X (jpl)                   
32 / 17 PC / 92 sc1bd4 [realization = 997946] / RUNNING (jpl)                  
33 / 17 PC / 93 sc1bd4 [realization = 991213] / C (jpl)                   
34 / 18 PC / 94 sc1bd4 [realization = 997946] / C (jpl)            
35 / 18 PC / 95 sc1bd4 [realization = 991213] / RUNNING (jpl)             
36 / 19 PC / 96 sc1bd4 [realization = 997946] / X (jpl)            
37 / 19 PC / 97 sc1bd4 [realization = 991213] / X (jpl)            
38 / 20 PC / 98 sc1bd4 [realization = 997946] / RUNNING (jpl)            
39 / 20 PC / 99 sc1bd4 [realization = 991213] / RUNNING (jpl)            

40 / 0 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)             
41 / 1 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)             
42 / 2 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)             
43 / 3 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)             
44 / 4 PC / sc1bd4 [realization = 9085] / waiting (jpl)            
45 / 5 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)             
46 / 6 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)             
47 / 7 PC / sc1bd4 [realization = 9085] / waiting (jpl)            
48 / 8 PC / sc1bd4 [realization = 9085] / waiting (jpl)            
49 / 9 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)             
50 / 10 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)             
51 / 11 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)             
52 / 12 PC / sc1bd4 [realization = 9085] / waiting (jpl)            
53 / 13 PC / sc1bd4 [realization = 9085] / waiting (jpl)            
54 / 14 PC / sc1bd4 [realization = 9085] / waiting (jpl)            
55 / 15 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)             
56 / 16 PC / sc1bd4 [realization = 9085] / waiting (jpl)            
57 / 17 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)             
58 / 18 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)             
59 / 19 PC / sc1bd4 [realization = 9085] / waiting (jpl)            
60 / 20 PC / sc1bd4 [realization = 9085] / RUNNING (jpl)                         

$\alpha\mathrm{s}\sim\mathrm{G}(\alpha_\mathrm{proj},0.003)$ [informative prior]  

OBS: The 0 PC is equivalent to MCMC40, so we can skip this MCMC.   
61 / 1 PC / sc1bd4 [realization = 9085] /  RUNNING (edge)                
62 / 2 PC / sc1bd4 [realization = 9085] /  RUNNING (edge)                
63 / 3 PC / sc1bd4 [realization = 9085] /  RUNNING (edge)                
64 / 4 PC / sc1bd4 [realization = 9085] /  RUNNING (edge)               
65 / 5 PC / sc1bd4 [realization = 9085] /  RUNNING (edge)                
66 / 6 PC / sc1bd4 [realization = 9085] /  RUNNING (edge)                
67 / 7 PC / sc1bd4 [realization = 9085] /  RUNNING (edge)               
68 / 8 PC / sc1bd4 [realization = 9085] /  RUNNING (edge)               
69 / 9 PC / sc1bd4 [realization = 9085] /  RUNNING (edge)                
70 / 10 PC / sc1bd4 [realization = 9085] / RUNNING (edge)                  
71 / 11 PC / sc1bd4 [realization = 9085] / RUNNING (edge)                  
72 / 12 PC / sc1bd4 [realization = 9085] / RUNNING (edge)                 
73 / 13 PC / sc1bd4 [realization = 9085] / RUNNING (edge)                 
74 / 14 PC / sc1bd4 [realization = 9085] / RUNNING (edge)                 
75 / 15 PC / sc1bd4 [realization = 9085] / RUNNING (edge)                  
76 / 16 PC / sc1bd4 [realization = 9085] / RUNNING (edge)                 
77 / 17 PC / sc1bd4 [realization = 9085] / RUNNING (edge)                  
78 / 18 PC / sc1bd4 [realization = 9085] / RUNNING (edge)                  
79 / 19 PC / sc1bd4 [realization = 9085] / RUNNING (edge)                 
80 / 20 PC / sc1bd4 [realization = 9085] / RUNNING (edge)       

!!! (So far, all these used cov_sc1bd4r0.) !!!   
!!! (From now, I'll use cov_sc1bd4 which was generated with nbar_sc1bd4) !!!     


$
\begin{align}
\boxed{\xi_\mathrm{fid}=\xi_\mathrm{cocoa}[n_{\mathrm{fid},\mathrm{DRM-D1}}^r]~;}
\end{align}
$     
= ======== Realizations under investigation ======== =  
Least extreme: r=3233  
Intermediate: r=0  
Most extreme: r=869  
= ======== Realizations under investigation ======== =  

Shift method: $\boxed{\xi_\mathrm{model}=\xi_\mathrm{cocoa}[\bar{n}_\mathrm{DRM-D1},\Delta_z^i]}$  
MCMC_BASELINEi / realization index / status       

22 / [realization = 3233] / converged (jpl)   
23 / [realization = 0] / converged (jpl)    
24 / [realization = 869] / converged (jpl)   
25 / [realization = 3233] / running (jpl)   
26 / [realization = 0] / running (jpl)    
27 / [realization = 869] / running (jpl)   
28 / [realization = 3233] / running (jpl)   
29 / [realization = 0] / running (jpl)    
30 / [realization = 869] / running (jpl)   

PCA method: $\boxed{\xi_\mathrm{model}=\xi_\mathrm{cocoa}[\bar{n}_\mathrm{DRM-D1}+\sum\alpha_i\mathrm{PC}_\mathrm{DRM-D1}]}$  
MCMCi / # of PCs / realization index / status   
$\alpha\mathrm{s}\sim\mathrm{G}(0,1)$ [uninformative prior]  

81 / 0 PC / [realization = 3233] / converged (jpl)    
82 / 1 PC / [realization = 3233] / converged (jpl)    
83 / 2 PC / [realization = 3233] / converged (jpl)    
84 / 3 PC / [realization = 3233] / converged (jpl)    
85 / 4 PC / [realization = 3233] / converged (jpl)    
86 / 5 PC / [realization = 3233] / converged (jpl)    
87 / 6 PC / [realization = 3233] /    
88 / 7 PC / [realization = 3233] /    
89 / 8 PC / [realization = 3233] /    
90 / 9 PC / [realization = 3233] /    
91 / 10 PC / [realization = 3233] /    
92 / 11 PC / [realization = 3233] /    
93 / 12 PC / [realization = 3233] /    
94 / 13 PC / [realization = 3233] /    
95 / 14 PC / [realization = 3233] /    
96 / 15 PC / [realization = 3233] /    
97 / 16 PC / [realization = 3233] /    
98 / 17 PC / [realization = 3233] /    
99 / 18 PC / [realization = 3233] /    
100 / 19 PC / [realization = 3233] /    
101 / 20 PC / [realization = 3233] /    

102 / 0 PC / [realization = 0] / converged (jpl)    
103 / 1 PC / [realization = 0] / converged (edge)    
104 / 2 PC / [realization = 0] / converged (edge)    
105 / 3 PC / [realization = 0] / converged (edge)    
106 / 4 PC / [realization = 0] / converged (edge)    
107 / 5 PC / [realization = 0] / converged (edge)    
108 / 6 PC / [realization = 0] /    
109 / 7 PC / [realization = 0] /    
110 / 8 PC / [realization = 0] /    
111 / 9 PC / [realization = 0] /    
112 / 10 PC / [realization = 0] /    
113 / 11 PC / [realization = 0] /    
114 / 12 PC / [realization = 0] /    
115 / 13 PC / [realization = 0] /    
116 / 14 PC / [realization = 0] /    
117 / 15 PC / [realization = 0] /    
118 / 16 PC / [realization = 0] /    
119 / 17 PC / [realization = 0] /    
120 / 18 PC / [realization = 0] /    
121 / 19 PC / [realization = 0] /    
122 / 20 PC / [realization = 0] /    

123 / 0 PC / [realization = 869] / converged (edge)    
124 / 1 PC / [realization = 869] / converged (edge)    
125 / 2 PC / [realization = 869] / converged (edge)    
126 / 3 PC / [realization = 869] / converged (edge)    
127 / 4 PC / [realization = 869] / converged (edge)    
128 / 5 PC / [realization = 869] / converged (edge)    
129 / 6 PC / [realization = 869] /    
130 / 7 PC / [realization = 869] /    
131 / 8 PC / [realization = 869] /    
132 / 9 PC / [realization = 869] /    
133 / 10 PC / [realization = 869] /    
134 / 11 PC / [realization = 869] /    
135 / 12 PC / [realization = 869] /    
136 / 13 PC / [realization = 869] /    
137 / 14 PC / [realization = 869] /    
138 / 15 PC / [realization = 869] /    
139 / 16 PC / [realization = 869] /    
140 / 17 PC / [realization = 869] /    
141 / 18 PC / [realization = 869] /    
142 / 19 PC / [realization = 869] /    
143 / 20 PC / [realization = 869] /    

$\alpha\mathrm{s}\sim\mathrm{G}(0,0.1)$ [slightly informative prior]  
$\alpha\mathrm{s}\sim\mathrm{G}(\alpha_\mathrm{proj},0.001)$ [informative prior]  


## PCs $\text{\textcolor{red}{trained}}$ on M1-D1 (sc1bd4), $\text{\textcolor{cyan}{reconstruct}}$ nz from M1-D2 (sc1bd5)   
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_ModelDRMD1_FidDRMD2/`  

$
\begin{align}
\boxed{\xi_\mathrm{model}=\xi_\mathrm{cocoa}[\bar{n}_\mathrm{DRM-D1}+\sum\alpha_i\mathrm{PC}_\mathrm{DRM-D1}]};
\end{align}
$

With shift $\Delta_{z,i}$. Sinthetic data: $\boxed{\xi_\mathrm{fid}=\xi_\mathrm{cocoa}[n_{\mathrm{fid},\mathrm{DRM-D2}}^r]~;~r=}$  

MCMC_BASELINEi / number of Δs per bin / status    
external_nz_modeling = 0  
1 / 1 /    
2 / 2 /      
3 / 3 /      
4 / 4 /      
5 / 5 /      
6 / 6 /      
7 / 7 /      
8 / 8 /      
9 / 9 /     

No shift ($\Delta_{z,i}=0$). Sinthetic data: $\boxed{\xi_\mathrm{fid}=\xi_\mathrm{cocoa}[n_{\mathrm{fid},\mathrm{DRM-D2}}^r]~;~r=}$  

0 / 0 PC / sc1bd5 [realization = ] / x (edge)          
1 / 1 PC / sc1bd5 [realization = ] / x (edge)          
2 / 2 PC / sc1bd5 [realization = ] / x (edge)          
3 / 3 PC / sc1bd5 [realization = ] / x (edge)          
4 / 4 PC / sc1bd5 [realization = ] / x (edge)          
5 / 5 PC / sc1bd5 [realization = ] / x (edge)          
6 / 6 PC / sc1bd5 [realization = ] / x (edge)          
7 / 7 PC / sc1bd5 [realization = ] / x (edge)          
8 / 8 PC / sc1bd5 [realization = ] / x (edge)          
9 / 9 PC / sc1bd5 [realization = ] / x (edge)          
10 / 10 PC / sc1bd5 [realization = ] / x (edge)                              
11 / 11 PC / sc1bd5 [realization = ] / x (edge)                    
12 / 12 PC / sc1bd5 [realization = ] / x (edge)                   
13 / 13 PC / sc1bd5 [realization = ] / x (edge)                   
14 / 14 PC / sc1bd5 [realization = ] / x (edge)                   
15 / 15 PC / sc1bd5 [realization = ] / x (edge)                    
16 / 16 PC / sc1bd5 [realization = ] / x (edge)                   
17 / 17 PC / sc1bd5 [realization = ] / x (edge)                    
18 / 18 PC / sc1bd5 [realization = ] / x (edge)                    
19 / 19 PC / sc1bd5 [realization = ] / x (edge)                   
20 / 20 PC / sc1bd5 [realization = ] / x (edge)              


## PCs $\text{\textcolor{red}{trained}}$ on M1-D2 (sc1bd5), $\text{\textcolor{cyan}{reconstruct}}$ nz from M1-D1 (sc1bd4)   
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_Model_sc1bd5_Fiducial_sc1bd4/`  

MCMCi  
0 / 0 PC / [realization = None] / converged (jpl)    
1 / 1 PC / [realization = None] / converged (jpl)    
2 / 2 PC / [realization = None] / converged (jpl)    
3 / 3 PC / [realization = None] / converged (jpl)    
4 / 4 PC / [realization = None] / converged (jpl)    
5 / 5 PC / [realization = None] / converged (jpl)    
MCMC_BASELINEi (0.003,0.01,0.1)      
0 / [realization = None] / running (jpl)       
1 / [realization = None] / running (jpl)             
2 / [realization = None] / running (jpl)        

## PCs $\text{\textcolor{red}{trained}}$ on M1-D3 (sc1bd6), $\text{\textcolor{cyan}{reconstruct}}$ nz from M1-D1 (sc1bd4)   
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_Model_sc1bd6_Fiducial_sc1bd4/`  

MCMCi  
0 / 0 PC / [realization = None] / converged (jpl)    
1 / 1 PC / [realization = None] / converged (jpl)    
2 / 2 PC / [realization = None] / converged (jpl)    
3 / 3 PC / [realization = None] / converged (jpl)    
4 / 4 PC / [realization = None] / converged (jpl)    
5 / 5 PC / [realization = None] / converged (jpl)    
MCMC_BASELINEi (0.003,0.01,0.1)     
0 / [realization = None] / running (jpl)          
1 / [realization = None] / running (jpl)            
2 / [realization = None] / running (jpl)            

## PCs $\text{\textcolor{red}{trained}}$ on M1-D4 (sc1bd7), $\text{\textcolor{cyan}{reconstruct}}$ nz from M1-D1 (sc1bd4)   
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_Model_sc1bd7_Fiducial_sc1bd4/`  

MCMCi  
0 / 0 PC / [realization = None] / converged (jpl)    
1 / 1 PC / [realization = None] / converged (jpl)    
2 / 2 PC / [realization = None] / converged (jpl)    
3 / 3 PC / [realization = None] / converged (jpl)    
4 / 4 PC / [realization = None] / converged (jpl)    
5 / 5 PC / [realization = None] / converged (jpl)     
MCMC_BASELINEi (0.003,0.01,0.1)       
0 / [realization = None] / running (jpl)         
1 / [realization = None] / running (jpl)            
2 / [realization = None] / running (jpl)            

## PCs $\text{\textcolor{red}{trained}}$ on M2-D2 (sc2bd5), $\text{\textcolor{cyan}{reconstruct}}$ nz from M2-D1 (sc2bd4)   
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_Model_sc2bd5_Fiducial_sc2bd4/`    
Transfer: `scp gattaca2-hn3:/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa/cocoa_photoz/results/chains/roman_pca/weighted_pca_Model_sc2bd5_Fiducial_sc2bd4/MCMC{0..5}* ./cocoa_photoz/results/chains/roman_pca/weighted_pca_Model_sc2bd5_Fiducial_sc2bd4/`   

MCMCi   
0 / 0 PC / [realization = None] / converged (edge)          
1 / 1 PC / [realization = None] / converged (edge)         
2 / 2 PC / [realization = None] / converged (edge)         
3 / 3 PC / [realization = None] / converged (edge)         
4 / 4 PC / [realization = None] / converged (edge)         
5 / 5 PC / [realization = None] / converged (edge)      
MCMC_BASELINEi (0.003,0.01,0.1)        
0 / [realization = None] / running (edge)            
1 / [realization = None] / running (edge)               
2 / [realization = None] / running (edge)               

## PCs $\text{\textcolor{red}{trained}}$ on M2-D3 (sc2bd6), $\text{\textcolor{cyan}{reconstruct}}$ nz from M2-D1 (sc2bd4)   
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_Model_sc2bd6_Fiducial_sc2bd4/`  
Transfer: `scp gattaca2-hn3:/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa/cocoa_photoz/results/chains/roman_pca/weighted_pca_Model_sc2bd6_Fiducial_sc2bd4/MCMC{0..5}* ./cocoa_photoz/results/chains/roman_pca/weighted_pca_Model_sc2bd6_Fiducial_sc2bd4/`   

MCMCi   
0 / 0 PC / [realization = None] / converged (edge)          
1 / 1 PC / [realization = None] / converged (edge)         
2 / 2 PC / [realization = None] / converged (edge)         
3 / 3 PC / [realization = None] / converged (edge)         
4 / 4 PC / [realization = None] / converged (edge)         
5 / 5 PC / [realization = None] / converged (edge)      
MCMC_BASELINEi (0.003,0.01,0.1)        
0 / [realization = None] / running (edge)  
1 / [realization = None] / running (edge)               
2 / [realization = None] / running (edge)               

## PCs $\text{\textcolor{red}{trained}}$ on M2-D4 (sc2bd7), $\text{\textcolor{cyan}{reconstruct}}$ nz from M2-D1 (sc2bd4)   
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_Model_sc2bd7_Fiducial_sc2bd4/`  
Transfer: `scp gattaca2-hn3:/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa/cocoa_photoz/results/chains/roman_pca/weighted_pca_Model_sc2bd7_Fiducial_sc2bd4/MCMC{0..5}* ./cocoa_photoz/results/chains/roman_pca/weighted_pca_Model_sc2bd7_Fiducial_sc2bd4/`   

MCMCi   
0 / 0 PC / [realization = None] / converged (edge)           
1 / 1 PC / [realization = None] / converged (edge)          
2 / 2 PC / [realization = None] / converged (edge)          
3 / 3 PC / [realization = None] / converged (edge)          
4 / 4 PC / [realization = None] / converged (edge)          
5 / 5 PC / [realization = None] / converged (edge)       
MCMC_BASELINEi (0.003,0.01,0.1)        
0 / [realization = None] / running (edge)          
1 / [realization = None] / running (edge)             
2 / [realization = None] / running (edge)             

## PCs $\text{\textcolor{red}{trained}}$ on M3-D2 (sc3bd5), $\text{\textcolor{cyan}{reconstruct}}$ nz from M3-D1 (sc3bd4)   
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_Model_sc3bd5_Fiducial_sc3bd4/`  

MCMCi   
0 / 0 PC / [realization = None] / converged (jpl)           
1 / 1 PC / [realization = None] / converged (jpl)          
2 / 2 PC / [realization = None] / converged (jpl)          
3 / 3 PC / [realization = None] / converged (jpl)          
4 / 4 PC / [realization = None] / converged (jpl)          
5 / 5 PC / [realization = None] / converged (jpl)       
MCMC_BASELINEi (0.003,0.01,0.1)        
0 / [realization = None] / running (jpl)            
1 / [realization = None] / running (edge)                     
2 / [realization = None] / running (edge)           

## PCs $\text{\textcolor{red}{trained}}$ on M3-D4 (sc3bd7), $\text{\textcolor{cyan}{reconstruct}}$ nz from M3-D1 (sc3bd4)   
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_Model_sc3bd7_Fiducial_sc3bd4/`  

MCMCi   
0 / 0 PC / [realization = None] / converged (jpl)          
1 / 1 PC / [realization = None] / converged (jpl)         
2 / 2 PC / [realization = None] / converged (jpl)         
3 / 3 PC / [realization = None] / converged (jpl)         
4 / 4 PC / [realization = None] / converged (jpl)         
5 / 5 PC / [realization = None] / converged (jpl)      
MCMC_BASELINEi (0.003,0.01,0.1)        
0 / [realization = None] / running (jpl)  
1 / [realization = None] / running (edge)                        
2 / [realization = None] / running (edge)              

# ###################################################################  
# ###################################################################  

# WEIGHTED PCA $w_0w_a\mathrm{CDM}$

## [DEPRECATED] PCs $\text{\textcolor{red}{trained}}$ on M1-D1 (sc1bd4), $\text{\textcolor{cyan}{reconstruct}}$ nz from M1-D1 (sc1bd4)    
(same as MCMC{61..80} of weighted PCA for $\Lambda\mathrm{CDM}$)  
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_Mod1d4_Fid1d4_w0wa/` - DEPRECATED  

$\alpha\mathrm{s}\sim\mathrm{G}(\alpha_\mathrm{proj},0.003)$ [informative prior]  

OBS: The 0 PC is equivalent to MCMC40, so we can skip this MCMC.   
1 / 1 PC / sc1bd4 [realization = 9085] /   paused (jpl)                 
2 / 2 PC / sc1bd4 [realization = 9085] /   paused (jpl)                 
3 / 3 PC / sc1bd4 [realization = 9085] /   paused (jpl)                 
4 / 4 PC / sc1bd4 [realization = 9085] /   paused (jpl)                
5 / 5 PC / sc1bd4 [realization = 9085] /   paused (jpl)                 
6 / 6 PC / sc1bd4 [realization = 9085] /   paused (jpl)                 
7 / 7 PC / sc1bd4 [realization = 9085] /   paused (jpl)                
8 / 8 PC / sc1bd4 [realization = 9085] /   paused (jpl)                
9 / 9 PC / sc1bd4 [realization = 9085] /   paused (jpl)                 
10 / 10 PC / sc1bd4 [realization = 9085] / paused (jpl)                   
11 / 11 PC / sc1bd4 [realization = 9085] / paused (jpl)                   
12 / 12 PC / sc1bd4 [realization = 9085] / paused (jpl)                  
13 / 13 PC / sc1bd4 [realization = 9085] / paused (jpl)                  
14 / 14 PC / sc1bd4 [realization = 9085] / paused (jpl)                  
15 / 15 PC / sc1bd4 [realization = 9085] / paused (jpl)                   
16 / 16 PC / sc1bd4 [realization = 9085] / paused (jpl)                  
17 / 17 PC / sc1bd4 [realization = 9085] / paused (jpl)                   
18 / 18 PC / sc1bd4 [realization = 9085] / paused (jpl)                   
19 / 19 PC / sc1bd4 [realization = 9085] / paused (jpl)                  
20 / 20 PC / sc1bd4 [realization = 9085] / paused (jpl)  

## PCs $\text{\textcolor{red}{trained}}$ on M1-D1 (sc1bd4), $\text{\textcolor{cyan}{reconstruct}}$ nz from M1-D1 (sc1bd4)   
(same as MCMC{81..86} of weighted PCA for $\Lambda\mathrm{CDM}$ but including +2 parameters: w0 and w0+wa)  
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_w0wa_Model_sc1bd4_Fiducial_sc1bd4/`    
Transfer: `scp gattaca2-hn3:/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa/cocoa_photoz/results/chains/roman_pca/weighted_pca_w0wa_Model_sc1bd4_Fiducial_sc1bd4/MCMC* ./cocoa_photoz/results/chains/roman_pca/weighted_pca_w0wa_Model_sc1bd4_Fiducial_sc1bd4/`   

MCMC_BASELINE      
0 / [realization = 3233] / x (jpl)        
1 / [realization = 0] / x (jpl)       
2 / [realization = 689] / x (jpl)       

MCMC  
0 / 0 PC / [realization = 3233] / running (edge)           
1 / 1 PC / [realization = 3233] / converged (edge)          
2 / 2 PC / [realization = 3233] / running (edge)          
3 / 3 PC / [realization = 3233] / converged (edge)          
4 / 4 PC / [realization = 3233] / running (edge)          
5 / 5 PC / [realization = 3233] / running (edge)       
...  
20 / 20 PC / [realization = 3233] / proxy

21 / 0 PC / [realization = 0] / running (edge)             
22 / 1 PC / [realization = 0] / running (edge)            
23 / 2 PC / [realization = 0] / converged (edge)            
24 / 3 PC / [realization = 0] / running (edge)            
25 / 4 PC / [realization = 0] / running (edge)            
26 / 5 PC / [realization = 0] / running (edge) - removed oversample and thin                         
...  
41 / 20 PC / [realization = 0] / proxy

42 / 0 PC / [realization = 869] / running (edge)            
43 / 1 PC / [realization = 869] / converged (edge)           
44 / 2 PC / [realization = 869] / running (edge) - removed oversample and thin                
45 / 3 PC / [realization = 869] / running (edge) - removed oversample and thin                
46 / 4 PC / [realization = 869] / running (edge) - removed oversample and thin                
47 / 5 PC / [realization = 869] / running (edge) - removed oversample and thin             
...  
62 / 20 PC / [realization = 869] / proxy


## PCs $\text{\textcolor{red}{trained}}$ on M1-D2 (sc1bd5), $\text{\textcolor{cyan}{reconstruct}}$ nz from M1-D1 (sc1bd4)   
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_w0wa_Model_sc1bd5_Fiducial_sc1bd4/`  

0 / 0 PC / [realization = None] / converged (jpl)    
1 / 1 PC / [realization = None] / x (jpl)    
2 / 2 PC / [realization = None] / running (jpl)    
3 / 3 PC / [realization = None] / x (jpl)    
4 / 4 PC / [realization = None] / converged (jpl)    
5 / 5 PC / [realization = None] / x (jpl) 

## PCs $\text{\textcolor{red}{trained}}$ on M1-D3 (sc1bd6), $\text{\textcolor{cyan}{reconstruct}}$ nz from M1-D1 (sc1bd4)   
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_w0wa_Model_sc1bd6_Fiducial_sc1bd4/`  

0 / 0 PC / [realization = None] / x (jpl)    
1 / 1 PC / [realization = None] / x (jpl)    
2 / 2 PC / [realization = None] / converged (jpl)    
3 / 3 PC / [realization = None] / x (jpl)    
4 / 4 PC / [realization = None] / x (jpl)    
5 / 5 PC / [realization = None] / x (jpl) 

## PCs $\text{\textcolor{red}{trained}}$ on M1-D4 (sc1bd7), $\text{\textcolor{cyan}{reconstruct}}$ nz from M1-D1 (sc1bd4)   
Yamls @ `cocoa_photoz/yamls/roman_pca/weighted_pca_w0wa_Model_sc1bd7_Fiducial_sc1bd4/`  

0 / 0 PC / [realization = None] / x (jpl)         
1 / 1 PC / [realization = None] / x (jpl)         
2 / 2 PC / [realization = None] / x (jpl)         
3 / 3 PC / [realization = None] / x (jpl)         
4 / 4 PC / [realization = None] / x (jpl)         
5 / 5 PC / [realization = None] / x (jpl)      