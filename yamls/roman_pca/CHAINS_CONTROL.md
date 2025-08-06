# cocoa_photoz  
MCMC PCA SCHEDULE PROGRAM   

# Model  
lcdm  
w0wa   

# Prior on alpha:  
Gaussian (0,1)  

# TEST 1 - fiducial cosmology from planck 18  
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

# TEST 2 - fiducial cosmology from planck 18  
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
36 sc1bd4 [realization index = 997946] / PAUSED (jpl)        [scp to edge done]  
37 sc1bd4 [realization index = 870731] / canceled (jpl)       [scp to edge done]    
38 sc1bd4 [realization index = 174820] / canceled (jpl)       [scp to edge done]    
39 sc1bd4 [realization index = 991213] / PAUSED (jpl)        [scp to edge done]    
40 sc1bd4 [realization index = 645431] / canceled (jpl)       [scp to edge done]    
41 sc1bd4 [realization index = 131103] / canceled (jpl)       [scp to edge done]    
42 sc1bd4 [realization index = 808898] / canceled (edge)      [scp to jpl done]      
43 sc1bd4 [realization index = 682038] / canceled (edge)      [scp to jpl done]      
44 sc1bd4 [realization index = 472580] / canceled (jpl)       [scp to edge done]    
45 sc1bd4 [realization index = 92550 ] / canceled (jpl)       [scp to edge done]     

### No PC + Shift with synthetic xi(nfid) + model xi(nbar)
46 sc1bd4 [realization index = 997946] / PAUSED (edge)       
47 sc1bd4 [realization index = 870731] / canceled (edge)      [scp to jpl done]   
48 sc1bd4 [realization index = 174820] / canceled (edge)      [scp to jpl done]     
49 sc1bd4 [realization index = 991213] / PAUSED (edge)         
50 sc1bd4 [realization index = 645431] / canceled (edge)      [scp to jpl done]     
51 sc1bd4 [realization index = 131103] / canceled (edge)      [scp to jpl done]     
52 sc1bd4 [realization index = 808898] / canceled (edge)      [scp to jpl done]     
53 sc1bd4 [realization index = 682038] / canceled (edge)      [scp to jpl done]     
54 sc1bd4 [realization index = 472580] / canceled (edge)      [scp to jpl done]     
55 sc1bd4 [realization index = 92550 ] / canceled (edge)      [scp to jpl done]     

### No PC + No Shift with synthetic xi(nfid) + model xi(fid)
56 sc1bd4 [realization index = 997946] / C (jpl)     
57 sc1bd4 [realization index = 991213] / C (jpl)     

### No PC + No Shift with synthetic xi(nfid) + model xi(nbar)
58 sc1bd4 [realization index = 997946] / C (jpl)     
59 sc1bd4 [realization index = 991213] / C (jpl)     

### 1 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
60 sc1bd4 [realization index = 997946] / C (jpl)     
61 sc1bd4 [realization index = 991213] / C (edge)     

### 2 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
62 sc1bd4 [realization index = 997946] / RUNNING (edge)     
63 sc1bd4 [realization index = 991213] / C (edge)     

### 3 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
64 sc1bd4 [realization index = 997946] / C (edge)     
65 sc1bd4 [realization index = 991213] / C (edge)     

### 4 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
66 sc1bd4 [realization index = 997946] / RUNNING (edge)     
67 sc1bd4 [realization index = 991213] / C (edge)     

### 5 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
68 sc1bd4 [realization index = 997946] / RUNNING (edge)     
69 sc1bd4 [realization index = 991213] / C (edge)     

### 6 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
70 sc1bd4 [realization index = 997946] / RUNNING (edge)     
71 sc1bd4 [realization index = 991213] / C (jpl)     

### 7 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
72 sc1bd4 [realization index = 997946] / RUNNING (jpl)     
73 sc1bd4 [realization index = 991213] / C (jpl)     

### 8 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
74 sc1bd4 [realization index = 997946] / WAITING (jpl)     
75 sc1bd4 [realization index = 991213] / C (jpl)     

### 9 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
76 sc1bd4 [realization index = 997946] / waiting (jpl)     
77 sc1bd4 [realization index = 991213] / RUNNING (jpl)     

### 10 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
78 sc1bd4 [realization index = 997946] / RUNNING (jpl)     
79 sc1bd4 [realization index = 991213] / C (jpl)       

### 11 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
80 sc1bd4 [realization index = 997946] / RUNNING (edge)     
81 sc1bd4 [realization index = 991213] / RUNNING (edge)      

### 12 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
82 sc1bd4 [realization index = 997946] / WAITING (edge)     
83 sc1bd4 [realization index = 991213] / WAITING (edge)      

### 13 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
84 sc1bd4 [realization index = 997946] / WAITING (edge)     
85 sc1bd4 [realization index = 991213] / WAITING (edge)      

### 14 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
86 sc1bd4 [realization index = 997946] / WAITING (edge)     
87 sc1bd4 [realization index = 991213] / WAITING (edge)      

### 15 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
88 sc1bd4 [realization index = 997946] / WAITING (edge)           
89 sc1bd4 [realization index = 991213] / WAITING (edge)      

### 16 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
90 sc1bd4 [realization index = 997946] / RUNNING (edge)      
91 sc1bd4 [realization index = 991213] / RUNNING (edge)      

### 17 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
92 sc1bd4 [realization index = 997946] / RUNNING (edge)     
93 sc1bd4 [realization index = 991213] / RUNNING (edge)      

### 18 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
94 sc1bd4 [realization index = 997946] / RUNNING (edge)     
95 sc1bd4 [realization index = 991213] / RUNNING (edge)      

### 19 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
96 sc1bd4 [realization index = 997946] / RUNNING (edge)     
97 sc1bd4 [realization index = 991213] / RUNNING (edge)      

### 20 PC + No Shift with synthetic xi(nfid) + model xi(nbar)
98 sc1bd4 [realization index = 997946] / RUNNING (edge)     
99 sc1bd4 [realization index = 991213] / RUNNING (edge)      

####################################################################  
####################################################################  



