Steps to run an evaluation and a mcmc with the pca model for n(z)
Assuming you're on ./cocoa/Cocoa

1) Copy and past the following files into the Roman folders: likelihood and interface
cp cocoa_photoz/codes/fisher.py projects/roman_real/likelihood/
cp cocoa_photoz/modifications/roman_real/_cosmolike_prototype_base.py projects/roman_real/likelihood/
cp cocoa_photoz/modifications/roman_real/roman_real_cosmic_shear.yaml projects/roman_real/likelihood/
cp cocoa_photoz/modifications/roman_real/roman_real_params_source.yaml projects/roman_real/likelihood/
cp cocoa_photoz/modifications/roman_real/interface.cpp projects/roman_real/interface/


2) Add the  dataset, lcdm.modelvector, mask, covariance, nz, U.txt ..? 

3) Compile Roman Real:
source projects/roman_real/scripts/compile_roman_real.sh 
