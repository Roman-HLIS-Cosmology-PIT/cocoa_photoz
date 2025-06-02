from cobaya.yaml import yaml_load_file
from cobaya.run import run
from cobaya.model import get_model


path = "/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz/"
info_from_yaml = yaml_load_file(path+"NZ_EVALUATE1.yaml")

info, _ = run(info_from_yaml)
model = get_model(info)
override_values = info["sampler"]["evaluate"]["override"]
model.logposterior(override_values)

likelihood = model.likelihood["lsst_y1.lsst_y1_cosmic_shear"]


print("likelihood.xi_diogo_test:",likelihood.dxi_dn)
