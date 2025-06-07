import sys, platform, os
os.environ['OMP_NUM_THREADS'] = '8'
import matplotlib
import math
from matplotlib import pyplot as plt
import numpy as np
import euclidemu2
import scipy

from getdist import IniFile
import itertools
import iminuit
import functools
print(sys.version)
print(os.getcwd())

# IMPORT CAMB
sys.path.insert(0, os.environ['ROOTDIR']+'/external_modules/code/CAMB/build/lib.linux-x86_64-'+os.environ['PYTHON_VERSION'])
import camb
from camb import model
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))

surveys = ['lsst_y1','roman_real']

for survey in surveys:

    if survey == 'lsst_y1':
        import cosmolike_lsst_y1_interface as ci
    elif survey == 'roman_real':    
        import cosmolike_roman_real_interface as ci

    print('\n---------------------------')
    print(f'cosmolike interface is {ci}')
    print('---------------------------\n')

    CAMBAccuracyBoost = 1.1
    non_linear_emul = 2
    CLprobe="xi"
    if survey == 'lsst_y1':
        path= "../external_modules/data/lsst_y1"
        data_file="lsst_y1_M1_GGL0.05.dataset"
    elif survey == 'roman_real':
        path= "../external_modules/data/roman_real"
        data_file="example1.dataset"

    IA_model = 0
    IA_redshift_evolution = 3

    ini = IniFile(os.path.normpath(os.path.join(path, data_file)))
    ci.initial_setup()
    ci.init_accuracy_boost(1.0, 1.0, int(1))
    ci.init_cosmo_runmode(is_linear = False)

    ci.init_redshift_distributions_from_files(
        lens_multihisto_file=ini.relativeFileName('nz_lens_file'),
        lens_ntomo=int(ini.int("lens_ntomo")), 
        source_multihisto_file=ini.relativeFileName('nz_source_file'),
        source_ntomo=int(ini.int("source_ntomo")))

    ci.init_IA( ia_model = int(IA_model), 
                ia_redshift_evolution = int(IA_redshift_evolution))

    ci.init_probes(possible_probes = CLprobe)

    ci.init_binning(int(ini.int("n_theta")), 
                    ini.float("theta_min_arcmin"), 
                    ini.float("theta_max_arcmin"))


    ci.init_data_real(ini.relativeFileName('cov_file'), 
                    ini.relativeFileName('mask_file'), 
                    ini.relativeFileName('data_file'))