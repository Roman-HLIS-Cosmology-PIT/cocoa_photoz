import os
os.environ['OMP_NUM_THREADS'] = '8'
import matplotlib.pyplot as plt
import numpy as np
import euclidemu2
from getdist import IniFile

# surveys = ['lsst_y1','roman_real','des_y3']
surveys = 'des_y3'


if survey == 'lsst_y1':
    import cosmolike_lsst_y1_interface as ci
elif survey == 'roman_real':    
    import cosmolike_roman_real_interface as ci
elif survey == 'des_y3':    
    import cosmolike_des_y3_interface as ci

print('\n---------------------------')
print(f'cosmolike interface is {ci}')
print('---------------------------\n')

CLprobe="xi"
if survey == 'lsst_y1':
    path= "../external_modules/data/lsst_y1"
    data_file="lsst_y1_M1_GGL0.05.dataset"
elif survey == 'des_y3':
    path= "../external_modules/data/des_y3"
    data_file="des_y3_real.dataset"
elif survey == 'roman_real':
    # path= "../external_modules/data/lsst_y1"
    # data_file="lsst_y1_M1_GGL0.05.dataset"
    path= "../external_modules/data/roman_real"
    data_file="example1.dataset"
    ggl_exclude = [[6,0],[7,0],[7,1]]
    ggl_exclude = np.array(ggl_exclude).flatten()

IA_model = 0
IA_redshift_evolution = 3

ini = IniFile(os.path.normpath(os.path.join(path, data_file)))
ci.initial_setup()
ci.init_accuracy_boost(1.0, 1.0, int(1))
ci.init_cosmo_runmode(is_linear = False)
if survey == 'roman_real':
    ci.init_ggl_exclude(ggl_exclude)
else: pass
ci.init_redshift_distributions_from_files(
    lens_multihisto_file=ini.relativeFileName('nz_lens_file'),
    lens_ntomo=int(ini.int("lens_ntomo")), 
    source_multihisto_file=ini.relativeFileName('nz_source_file'),
    source_ntomo=int(ini.int("source_ntomo")))

ci.init_IA( ia_model = int(IA_model), 
            ia_redshift_evolution = int(IA_redshift_evolution))

ci.init_probes(possible_probes = CLprobe)

x = ci.init_binning(int(ini.int("n_theta")), 
                ini.float("theta_min_arcmin"), 
                ini.float("theta_max_arcmin"))

ci.init_data_real(ini.relativeFileName('cov_file'), 
                ini.relativeFileName('mask_file'), 
                ini.relativeFileName('data_file'))

x = np.array(ci.get_dv_masked())
# np.savetxt("./get_dv_masked.txt",x)
